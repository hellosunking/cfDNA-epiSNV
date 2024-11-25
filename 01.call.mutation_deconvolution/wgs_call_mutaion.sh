#!/bin/bash

sample_list=$1

# Check if the sample list file exists
if [[ ! -f $sample_list ]]; then
  echo "Error: Sample list file '$sample_list' not found!"
  exit 1
fi

# Iterate over each sample in the list
while read -r sid group; do
  echo "Processing sample: $sid"
  
  # Step 1: Filter BAM file
  samtools view -h $sid/$sid.mkdup.bam | \
  perl -alne 'print and next if /^@/; next if $F[1] & 0x100; next unless $F[1] & 0x02; next unless $F[6] eq "="; next if $F[2] eq "*"; print' | \
  samtools view -b -o $sid/$sid.filter.bam

  samtools index $sid/$sid.filter.bam

  # Step 2: Generate VCF file
  bcftools mpileup -q 60 -Q 30 -Ov -f /usr/genomes/hg38p13.fa $sid/$sid.filter.bam | \
  bcftools call -vc -V indels -Oz -o $sid/$sid.vcf.gz

  tabix $sid/$sid.vcf.gz

  # Step 3: Filter VCF
  zcat $sid/$sid.vcf.gz | \
  perl -alne 'print and next if /^#/; print if $F[0] ne "chrX" && $F[0] ne "chrY" && $F[0] ne "chrM" && $F[5] >= 30.0 && $F[0] !~ /_/ && $F[4] !~ /,/' | \
  bgzip > $sid/$sid.hdfilt.pass.snp.raw.vcf.gz

  zcat $sid/$sid.hdfilt.pass.snp.raw.vcf.gz | \
  perl -alne 'next if /^#/; print join "\t", $F[0], $F[1]' > $sid/$sid.pos1.txt

  perl get_dmbs_pos.pl $sid/$sid.pos1.txt | uniq > $sid/$sid.dmbs1.txt

  perl filt_vcf_dmbs.pl $sid/$sid.hdfilt.pass.snp.raw.vcf.gz $sid/$sid.dmbs1.txt $sid/$sid.hdfilt.pass.snp.tmp.vcf.gz

  tabix $sid/$sid.hdfilt.pass.snp.tmp.vcf.gz

  # Step 4: Exclude dark regions
  bcftools view $sid/$sid.hdfilt.pass.snp.tmp.vcf.gz -R /usr/genomes/hg38.remove.Darkregion.bed | \
  bgzip > $sid/$sid.hdfilt.pass.snp.tmp1.vcf.gz

  perl sub_dp.pl $sid/$sid.hdfilt.pass.snp.tmp1.vcf.gz $sid/$sid.hdfilt.pass.snp.tmp2.vcf.gz

  # Step 5: Final filtering
  perl filt_vcf_snp.pl $sid/$sid.hdfilt.pass.snp.tmp2.vcf.gz /usr/genomes/hg38_cosmic99.filt.txt /usr/genomes/known_snp.txt $sid/$sid.final.vcf.gz

  tabix $sid/$sid.final.vcf.gz

  echo "Sample $sid processing complete."
done < "$sample_list"

