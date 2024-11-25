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
  
  # Step 1: Deal with BAM file
  samtools view -@ 2 -h $sid/Msuite2.final.bam | perl -alne 'print and next if /^@/; next if $F[1] & 0x100 ; next unless $F[1] & 0x02; next unless $F[6] eq "="; next if $F[2] eq "*"; print if /XG:Z:CT/' | samtools view -@ 2 -b -o $sid/$sid.ct.bam &

  samtools view -@ 2 -h $sid/Msuite2.final.bam | perl -alne 'print and next if /^@/; next if $F[1] & 0x100 ; next unless $F[1] & 0x02; next unless $F[6] eq "="; next if $F[2] eq "*"; print if /XG:Z:GA/' | samtools view -@ 2 -b -o $sid/$sid.ga.bam &
  wait

  samtools index -@ 2 $sid/$sid.ct.bam &
  samtools index -@ 2 $sid/$sid.ga.bam &
  wait

  # Step 2: Generate VCF file
  bcftools mpileup --threads 2 -q 60 -Q 30 -Ov -f /usr/genomes/hg38_p13/hg38p13.fa $sid/$sid.ct.bam | bcftools call --threads 2 -vc -V indels -Oz -o $sid/c.vcf.gz &

  bcftools mpileup --threads 2 -q 60 -Q 30 -Ov -f /usr/genomes/hg38_p13/hg38p13.fa $sid/$sid.ga.bam | bcftools call --threads 2 -vc -V indels -Oz -o $sid/g.vcf.gz &
  wait

  zcat $sid/c.vcf.gz | perl -alne 'print and next if /^#/; print if $F[4] !~ /T/' > $sid/c1.vcf &
  zcat $sid/g.vcf.gz | perl -alne 'print and next if /^#/; print if $F[4] !~ /A/' > $sid/g1.vcf &
  wait

  bgzip ../02.filter/$sid/c1.vcf &
  bgzip ../02.filter/$sid/g1.vcf &
  wait

  tabix $sid/c1.vcf.gz &
  tabix $sid/g1.vcf.gz &
  wait

  bcftools merge --threads 4 $sid/c1.vcf.gz $sid/g1.vcf.gz -Oz -o $sid/$sid.vcf.gz
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

