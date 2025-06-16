## Codes and scripts used in Zhang and An et al. manuscript
Distributed under the [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/ "CC BY-NC-ND") license and for **personal and academic usage only**.

The following softwares are used:
- [ktrim v1.5.0](https://github.com/hellosunking/Ktrim/releases/tag/v1.5.0)
- [bwa v0.7.17](https://github.com/lh3/bwa/releases/tag/v0.7.17)
- [GATK v4.2.3.0](https://github.com/broadinstitute/gatk/releases/tag/4.2.3.0)
- [picard v2.23.4](https://github.com/broadinstitute/picard/releases/tag/2.23.4)
- [Msuite2 v2.1.0](https://github.com/hellosunking/Msuite2/releases/tag/v2.1.0)
- [samtools v1.17](https://github.com/samtools/samtools/releases/tag/1.17)
- [bcftools v1.12](https://github.com/samtools/bcftools/releases/tag/1.12)
- [R v4.2.0](https://cran.r-project.org/bin/windows/base/old/4.2.0), requires "ggplot2", "ggpubr", "ComplexHeatmap", "MutationalPatterns", "BSgenome", "NMF", "gbm", "caret", "foreach", "doParallel", "binom", and "pROC" packages.
- [bgzip v1.17](http://www.htslib.org/doc/1.17/bgzip.html)
- [tabix v1.17](http://www.htslib.org/doc/1.17/tabix.html)
- [bedtools v2.29.2](https://github.com/arq5x/bedtools2/releases/tag/v2.29.2)

To run the analysis, you need to clone this repo to obtain all the programs and related files:
```
git clone https://github.com/hellosunking/cfDNA-epiSNV.git
## you will see a newly created directory named "cfDNA-epiSNV", go into it for following analysis
cd cfDNA-epiSNV
```

We had deposited processed data (including variants, cfDNA reads, and large annotation files) to [Zenodo](https://zenodo.org/records/14849892 "Zenodo data link"). You may download it and decompress the files at the same directory of this README file:
```
wget -O fragmentomics.of.variants.in.cfDNA.tar "https://zenodo.org/records/14849892/files/fragmentomics.of.variants.in.cfDNA.tar?download=1"
tar xf fragmentomics.of.variants.in.cfDNA.tar
wget -O suppl.tar "https://zenodo.org/records/XXXXXXXXXX/files/suppl.CH.tar?download=1"
tar xf suppl.CH.tar
## you will see a newly created directory named "Processed.files" with all files stored inside.
```

## 1. Read preprocessing and alignment
```
## Need to update sampleID and path to the FASTQ files
sid=SampleID
FASTQ1=/path/to/$sid.R1.fq.gz
FASTQ2=/path/to/$sid.R2.fq.gz

#### For plain DNA-seq data
## data preprocess. Note that all the WGS data we analyzed are generated using MGI sequencers
ktrim -1 $FASTQ1 -2 $FASTQ2 -t 8 -o $sid.ktrim -k BGI

## read alignment. Need to build hg38 index for BWA first and replace the following path
hg38index=/path/to/hg38.bwa.index
bwa mem -t 16 -M -Y $hg38index $sid.ktrim.read1.fq $sid.ktrim.read2.fq | samtools view -b -@ 16 - | samtools sort -@ 16 -o $sid.sort.bam - 2> $sid.bwa.log

## remove duplicates
java -Xmx16g -XX:ParallelGCThreads=8 -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$sid.sort.bam O=$sid.mkdup.bam METRICS_FILE=$sid.mkdup.bam.mat TMP_DIR=tmp 2> $sid.mkdup.log

#### For EM-seq data
## Run Msuite2. Need to build hg38 index for Msuite2 first. The data is sequenced on illumina sequencers
msuite2 -x hg38 -1 $FASTQ1 -2 $FASTQ2 -o Msuite2.$sid -k illumina --cut-r1-tail 25 --cut-r2-head 25 --aligner hisat2 -p 16
cd Msuite2.$sid
make
make clean
cd ../
```

## 2. Identification of CH- and tumor-derived somatic variants in paired PBMC and tumor genotyping data
```
## preprocess data, read alignment and BQSR are the same as Step 1
sid=BRCA
## assumes that you have already get  $sid.PBMC.mkdup.bam and $sid.tumor.mkdup.bam

## BQSR using GATK
## Need to update the path to hg38 genome file (FASTA format) 
hg38fasta=/path/to/hg38.fa

gatk BaseRecalibrator -R $hg38fasta -I $sid.PBMC.mkdup.bam -O $sid.PBMC.recal.table \
	--known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
	--known-sites 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
	--known-sites dbsnp_156.vcf.gz
gatk ApplyBQSR -R $hg38fasta -I $sid.PBMC.mkdup.bam -O $sid.PBMC.recal.bam \
	-bqsr-recal-file $sid.PBMC.recal.table

gatk BaseRecalibrator -R $hg38fasta -I $sid.tumor.mkdup.bam -O $sid.tumor.recal.table \
	--known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
	--known-sites 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
	--known-sites dbsnp_156.vcf.gz
gatk ApplyBQSR -R $hg38fasta -I $sid.tumor.mkdup.bam -O $sid.tumor.recal.bam \
	-bqsr-recal-file $sid.tumor.recal.table

## call variants using GATK mutect2
gatk Mutect2 -R $hg38fasta --native-pair-hmm-threads 32 \
		-I $sid.tumor.recal.bam \
		-I $sid.PBMC.recal.bam \
		-L Agilent.SureSelectXT.Human.All.Exon.V4+UTRs.bed \
		-normal $sid.PBMC \
		--genotype-germline-sites \
		--min-base-quality-score 20 \
		--read-filter MappingQualityReadFilter \
		--minimum-mapping-quality 20 \
		--dont-use-soft-clipped-bases \
		-O $sid.tumor.vcf.gz
gatk SelectVariants -V $sid.PBMC.vcf.gz  -select-type SNP -O $sid.PBMC.SNP.vcf.gz

gatk Mutect2 -R $hg38fasta --native-pair-hmm-threads 32 \
		-I $sid.tumor.recal.bam \
		-I $sid.PBMC.recal.bam \
		-L Agilent.SureSelectXT.Human.All.Exon.V4+UTRs.bed \
		-normal $sid.tumor \
		--genotype-germline-sites \
		--min-base-quality-score 20 \
		--read-filter MappingQualityReadFilter \
		--minimum-mapping-quality 20 \
		--dont-use-soft-clipped-bases \
		-O $sid.PBMC.vcf.gz
gatk SelectVariants -V $sid.tumor.vcf.gz -select-type SNP -O $sid.tumor.SNP.vcf.gz

## filter variants, PBMC-specific variants are considered CH-derived
## the processed files in vcf format are under Processed.files/6.Exome-seq/ directory.
perl 1.variant/filt_paired_pbmc.pl  $sid.PBMC.vcf.gz  | bgzip > $sid.PBMC.filter.vcf.gz
perl 1.variant/filt_paired_tumor.pl $sid.tumor.vcf.gz | bgzip > $sid.tumor.filter.vcf.gz
perl 1.variant/filt_paired_vcf.pl $sid.PBMC.filter.vcf.gz $sid.tumor.filter.vcf.gz | bgzip > $sid.CH.vcf.gz
perl 1.variant/filt_paired_vcf.pl $sid.tumor.filter.vcf.gz $sid.PBMC.filter.vcf.gz | bgzip > $sid.tumor.vcf.gz
## the processed files in vcf format are under Processed.files/6.CH.vs.Tumor/ directory.
```

## 3. Identification of CH-dervied variants in high-depth PBMC data
```
## here we use the sample Ctrl_1 as an example
## You need to uncompress the SRA file use fasterq-dump to get fastq files
sid=Ctrl_1
PBMC_1=SRR10799887_1.fastq
PBMC_2=SRR10799887_2.fastq

## extract and remove UMIs in this data, using the processed FASTQ files for read alignment
perl 4.target-seq/process.UMI.pl $PBMC_1 $PBMC_2 $sid.PBMC

## preprocessing, read alignment, and remove duplicates, same as Step 1
## and assumes that you get the file $sid.PBMC.mkdup.bam
## call variants using GATK, the annotation files for GATK are required
## Need to update the path to hg38 genome file (FASTA format) 
hg38fasta=/path/to/hg38.fa

gatk BaseRecalibrator -R $hg38fasta -I $sid.PBMC.mkdup.bam -O $sid.PBMC.recal.table \
	--known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
	--known-sites 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
	--known-sites dbsnp_156.vcf.gz

gatk ApplyBQSR -R $hg38fasta -I $sid.PBMC.mkdup.bam -O $sid.PBMC.recal.bam \
	-bqsr-recal-file $sid.PBMC.recal.table

gatk HaplotypeCaller -R $hg38fasta -I $sid.PBMC.recal.bam -O $sid.PBMC.vcf.gz \
	-ERC GVCF --dont-use-soft-clipped-bases true --min-base-quality-score 20 \
	--read-filter MappingQualityReadFilter --minimum-mapping-quality 20

gatk GenotypeGVCFs -R $hg38fasta -V $sid.PBMC.vcf.gz \
	-D dbsnp_156.vcf.gz -O $sid.PBMC.HC.vcf.gz

gatk SelectVariants -V $sid.PBMC.HC.vcf.gz -select-type SNP -O $sid.PBMC.SNP.vcf.gz

## filter variants
## the processed files in vcf format are under Processed.files/5.Target-seq/ directory.
perl 4.target-seq/filt_vcf_pbmc.pl $sid.PBMC.SNP.vcf.gz $sid.CH.vcf.gz
```

## 4. Call and filter variants in low-pass cfDNA data
```
#### For plain DNA-seq data:
## Need to update the path to hg38 genome file (FASTA format)
hg38fasta=/path/to/hg38.fa
samtools view -@ 4 -h $sid.mkdup.bam | perl -alne 'print and next if /^@/; next if $F[1] & 0x100 || $F[2] eq "*"; print if ($F[1] & 0x02) && ($F[6] eq "=")' | samtools view -b -@ 16 -o $sid.filter.bam
samtools index -@ 16 $sid.filter.bam
bcftools mpileup -q 60 -Q 30 -Ov -f $hg38fasta --threads 16 $sid.filter.bam | bcftools call -vc -V indels -Oz --threads 16 -o $sid.vcf.gz
zcat $sid.vcf.gz | perl -alne 'print and next if /^#/; print if $F[0]=~/^chr\d+$/ && $F[5] >= 30 && $F[4] !~ /,/' | gzip > $sid.hdfilt.pass.snp.raw.vcf.gz
zcat $sid.hdfilt.pass.snp.raw.vcf.gz | perl -alne 'next if /^#/; print "$F[0]\t$F[1]"' | perl 1.variant/get_dmbs_pos.pl - | uniq > $sid.dmbs.txt
perl 1.variant/filt_vcf_dmbs.pl $sid.hdfilt.pass.snp.raw.vcf.gz $sid.dmbs.txt $sid.hdfilt.pass.snp.tmp.vcf.gz

## the hg38.non.blacklist.bed records the genomic regions that are NOT in ENCODE's blacklist
bcftools view $sid.hdfilt.pass.snp.tmp.vcf.gz -R 1.variant/hg38.non.blacklist.bed | gzip > $sid.hdfilt.pass.snp.tmp1.vcf.gz
perl 1.variant/sub_dp.pl $sid.hdfilt.pass.snp.tmp1.vcf.gz $sample.hdfilt.pass.snp.tmp2.vcf.gz

## the "hg38_cosmic99.filt.gz" and "known_snp.gz" files could be found in "Processed.files"
perl 1.variant/filt_vcf_snp.pl $sid.hdfilt.pass.snp.tmp2.vcf.gz Processed.files/hg38_cosmic99.filt.gz Processed.files/known_snp.gz $sid.final.vcf.gz

#### For EM-seq data
samtools view -@ 4 -h Msuite2.$sid/Msuite2.final.bam | perl -alne 'print and next if /^@/; next if $F[1] & 0x100 || $F[2] eq "*"; print if ($F[1] & 0x02) && ($F[6] eq "=") && /XG:Z:CT/' | samtools view -@ 8 -b -o $sid.ct.bam &
samtools view -@ 4 -h Msuite2.$sid/Msuite2.final.bam | perl -alne 'print and next if /^@/; next if $F[1] & 0x100 || $F[2] eq "*"; print if ($F[1] & 0x02) && ($F[6] eq "=") && /XG:Z:GA/' | samtools view -@ 8 -b -o $sid.ga.bam &
wait

bcftools mpileup --threads 4 -q 60 -Q 30 -Ov -f $hg38fasta $sid.ct.bam | bcftools call --threads 4 -vc -V indels -Oz -o $sid.c.vcf.gz &
bcftools mpileup --threads 4 -q 60 -Q 30 -Ov -f $hg38fasta $sid.ga.bam | bcftools call --threads 4 -vc -V indels -Oz -o $sid.g.vcf.gz &
wait

## filter T variants, specific to EM-seq data
zcat $sid.c.vcf.gz | perl -alne 'print and next if /^#/; print if $F[4] !~ /T/' | bgzip > $sid.c.filter.vcf.gz &
zcat $sid.g.vcf.gz | perl -alne 'print and next if /^#/; print if $F[4] !~ /A/' | bgzip > $sid.g.filter.vcf.gz &
wait
tabix $sid.c.filter.vcf.gz &
tabix $sid.g.filter.vcf.gz &
wait
bcftools merge --threads 16 $sid.c.filter.vcf.gz $sid.g.filter.vcf.gz -Oz -o $sid.vcf.gz
```

For each dataset, after obtaining the somatic variants for all samples, the following commands were used to generate mutation profiles:
```
## prepare a "vcf.info" file, which contains 3 columns: sampleID category /path/to/final.vcf.gz
## we provided example "vcf.info" files for HCC, Liang, and Bie cohorts
## The following is the codes to analyze HCC cohort samples
cohortID=HCC
Rscript 1.variant/MutationalPatterns_profile.R 1.variant/$cohortID.vcf.info 1.variant/$cohortID.final.mut_mat.txt
perl 1.variant/normalize_matrix.pl 1.variant/$cohortID.final.mut_mat.txt | awk '{for(i=1;i<=NF;i++)a[i]=a[i]?a[i]"\t"$i:$i}END{for(i=1;i<=NF;i++)print a[i]}' > 1.variant/$cohortID.96_mutation_profile.txt
## the 96_mutation_profile.txt files for all cohorts are also included in this package

## PCA and cluster for HCC cohort
Rscript 1.variant/pca.R 1.variant/$cohortID.96_mutation_profile.txt 1.variant/$cohortID.vcf.info $cohortID
Rscript 1.variant/cluster.R 1.variant/$cohortID.96_mutation_profile.txt 1.variant/$cohortID.vcf.info $cohortID
```

## 5. Extract Wt- and Mut-DNA
For each sample, we used the following commands to extract Wt- and Mut- DNA in BED format (which files are provided in "Processed.files").
```
zcat $sid.final.vcf.gz | perl -alne 'next if /^#/; print join "\t", $F[0], $F[1], $F[3], $F[4], $F[2]' > $sid.txt

## set EMseq=0 for WGS data, and EMseq=1 for EM-seq data
EMseq=0

while read chr locus ref mut extra
do
    ## extend the loci to grab reads
    let left=$locus-1000
    let right=$locus+1000

    if [ $EMseq == "0" ] ## WGS-seq data
    then
        samtools view $sid.mkdup.bam "$chr:$left-$right" | perl 2.fragmentomics_epigenetics/1.extract_fragments/process.sam_wgs.pl - "$chr:$locus" "$locus,$ref,$mut" 1
    else    ## EM-seq data
        samtools view $sid.mkdup.bam "$chr:$left-$right" | perl 2.fragmentomics_epigenetics/1.extract_fragments/process.sam_EM.pl  - "$chr:$locus" "$locus,$ref,$mut" 1
    fi
done < $sid.txt > $sid.sequenced.alleles.tmp

perl 2.fragmentomics_epigenetics/1.extract_fragments/select.pl $sid.sequenced.alleles.tmp | perl 2.fragmentomics_epigenetics/1.extract_fragments/rm_more_loci.pl - > $sid.sequenced.alleles

## get the Mut- and Wt-DNA files in BED format
sort -k1,1V -k2,2n $sid.sequenced.alleles | perl 2.fragmentomics_epigenetics/1.extract_fragments/extract.bed.pl - $sid
```

## 6. CfDNA fragmentomic features in Wt- and Mut-DNA
After extracting Wt- and Mut-DNA, we used the following commands to analyze cfDNA fragmentomic features for each sample:
```
## All processed Wt- and Mut-DNA files are under "Processed.files"
## Need to update the file paths
WtDNA=Processed.files/1.HCC/Wt_DNA/$sid.bed.gz
MutDNA=Processed.files/1.HCC/Mut_DNA/$sid.bed.gz

## size
perl 2.fragmentomics_epigenetics/2.size/bed2size.pl $WtDNA  $sid.Wt.size
perl 2.fragmentomics_epigenetics/2.size/bed2size.pl $MutDNA $sid.Mut.size

Rscript 2.fragmentomics_epigenetics/2.size/plot_size.R $sid
n=`perl 2.fragmentomics_epigenetics/2.size/size_freq.pl $sid.Wt.size`
t=`perl 2.fragmentomics_epigenetics/2.size/size_freq.pl $sid.Mut.size`
echo -e "$sid\t$n\t$t" | perl -alne 'print join "\t", $F[0], $F[1], $F[2], $F[2]-$F[1]' > $sid.Diff_size

## end motif. Need to update the path to hg38 genome file (FASTA format)
hg38fasta=/path/to/hg38.fa
perl 2.fragmentomics_epigenetics/3.motif/grab.end.with.extension.pl $WtDNA  $sid.Wt  PE &
perl 2.fragmentomics_epigenetics/3.motif/grab.end.with.extension.pl $MutDNA $sid.Mut PE &
wait
perl 2.fragmentomics_epigenetics/3.motif/bed2fa.local.pl $hg38fasta $sid.Wt.left.outer2.inner4.bed $sid.Mut.left.outer2.inner4.bed
perl 2.fragmentomics_epigenetics/3.motif/extract.motif.with.extension.pl $sid.Wt.left.outer2.inner4.fa  > $sid.Wt.motif  &
perl 2.fragmentomics_epigenetics/3.motif/extract.motif.with.extension.pl $sid.Mut.left.outer2.inner4.fa > $sid.Mut.motif &
wait
## clean the intermediate files (optional)
#rm -f $sid.Wt.left.outer2.inner4.bed $sid.Mut.left.outer2.inner4.bed $sid.Wt.left.outer2.inner4.fa $sid.Mut.left.outer2.inner4.fa

ccca_n=$(grep -w "^CCCA" $sid.Wt.motif  | cut -f 5)
ccca_t=$(grep -w "^CCCA" $sid.Mut.motif | cut -f 5)
echo -e "$sid\t$ccca_n\t$ccca_t" | perl -alne 'print join "\t", $F[0], $F[1], $F[2], $F[2]-$F[1]' > $sid.motif.ccca

ctcc_n=$(grep -w "^CTCC" $sid.Wt.motif  | cut -f 3)
ctcc_t=$(grep -w "^CTCC" $sid.Mut.motif | cut -f 3)
echo -e "$sid\t$ctcc_n\t$ctcc_t" | perl -alne 'print join "\t", $F[0], $F[1], $F[2], $F[2]-$F[1]' > $sid.motif.ctcc

## nucleosome_track
bedtools intersect -a $WtDNA  -b Processed.files/GM12878.nucleosome.hg38.bed.gz -sorted -wao | perl 2.fragmentomics_epigenetics/4.nucleosome/anno.end.pl $sid.Wt  - &
bedtools intersect -a $MutDNA -b Processed.files/GM12878.nucleosome.hg38.bed.gz -sorted -wao | perl 2.fragmentomics_epigenetics/4.nucleosome/anno.end.pl $sid.Mut - &
wait

Rscript 2.fragmentomics_epigenetics/4.nucleosome/plot.end.R $sid
perl 2.fragmentomics_epigenetics/4.nucleosome/calc_nuc50.pl $sid > $sid.Diff_nuc

## E-index. Need to update the path to the end model, which is extremely large (~12GB).
## You can build it following the instructions in https://github.com/hellosunking/molecular-cfDNA-fragmentomics,
## or approach the corresponding author for the file.
EindexEndModel=/path/to/EindexEndModel
echo -e "$sid.Wt\t$WtDNA"    > $sid.bed.list
echo -e "$sid.Mut\t$MutDNA" >> $sid.bed.list
2.fragmentomics_epigenetics/5.Eindex/calc.E-index.multi-thread 2.fragmentomics_epigenetics/5.Eindex/hg38.genome.info $EindexEndModel 2.fragmentomics_epigenetics/5.Eindex/ENCODE.blacklist.hg38.bed $sid.bed.list 2 y > $sid.Eindex.raw
perl 2.fragmentomics_epigenetics/5.Eindex/calc.Diff.Eindex.pl $sid > $sid.Diff.Eindex
```

For each cohort, after calculating the fragmentomics features for all samples, we sorted the results and used the following commands to do visualizations:
```
## use HCC cohort as an example. The sorted results for this cohort are provided in this package
cohortID=HCC
## CancerType could be ESCA, COREAD, etc for Bie et al. cohort
CancerType=HCC
ControlType=Control

## size
Rscript 2.fragmentomics_epigenetics/line.R 2.fragmentomics_epigenetics/2.size/$cohortID.Diff_size $cohortID.$CancerType.S150 $CancerType 2.fragmentomics_epigenetics/2.size/$cohortID.$CancerType.line.pdf
Rscript 2.fragmentomics_epigenetics/line.R 2.fragmentomics_epigenetics/2.size/$cohortID.Diff_size $cohortID.$ControlType.S150 $ControlType 2.fragmentomics_epigenetics/2.size/$cohortID.$ControlType.line.pdf

Rscript 2.fragmentomics_epigenetics/box.R 2.fragmentomics_epigenetics/2.size/$cohortID.Diff_size Diff_size 2.fragmentomics_epigenetics/2.size/$cohortID.Diff_size.box.pdf
Rscript 2.fragmentomics_epigenetics/roc.R 2.fragmentomics_epigenetics/2.size/$cohortID.Diff_size Diff_size 2.fragmentomics_epigenetics/2.size/$cohortID.Diff_size.roc.pdf

## end motif
Rscript 2.fragmentomics_epigenetics/line.R 2.fragmentomics_epigenetics/3.motif/$cohortID.Diff_CCCA $CancerType.CCCA $CancerType 2.fragmentomics_epigenetics/3.motif/$cohortID.$CancerType.ccca_line.pdf
Rscript 2.fragmentomics_epigenetics/line.R 2.fragmentomics_epigenetics/3.motif/$cohortID.Diff_CCCA $ControlType.CCCA $ControlType 2.fragmentomics_epigenetics/3.motif/$cohortID.$ControlType.ccca_line.pdf
Rscript 2.fragmentomics_epigenetics/box.R 2.fragmentomics_epigenetics/3.motif/$cohortID.Diff_CCCA Diff_CCCA 2.fragmentomics_epigenetics/3.motif/$cohortID.Diff_CCCA.box.pdf
Rscript 2.fragmentomics_epigenetics/roc.R 2.fragmentomics_epigenetics/3.motif/$cohortID.Diff_CCCA Diff_CCCA 2.fragmentomics_epigenetics/3.motif/$cohortID.Diff_CCCA.roc.pdf

Rscript 2.fragmentomics_epigenetics/line.R 2.fragmentomics_epigenetics/3.motif/$cohortID.Diff_CTCC $CancerType.CTCC $CancerType 2.fragmentomics_epigenetics/3.motif/$cohortID.$CancerType.ctcc_line.pdf
Rscript 2.fragmentomics_epigenetics/line.R 2.fragmentomics_epigenetics/3.motif/$cohortID.Diff_CTCC $ControlType.CTCC $ControlType 2.fragmentomics_epigenetics/3.motif/$cohortID.$ControlType.ctcc_line.pdf
Rscript 2.fragmentomics_epigenetics/box.R 2.fragmentomics_epigenetics/3.motif/$cohortID.Diff_CTCC Diff_CTCC 2.fragmentomics_epigenetics/3.motif/$cohortID.Diff_CTCC.box.pdf
Rscript 2.fragmentomics_epigenetics/roc.R 2.fragmentomics_epigenetics/3.motif/$cohortID.Diff_CTCC Diff_CTCC 2.fragmentomics_epigenetics/3.motif/$cohortID.Diff_CTCC.roc.pdf

## nucleosome track
Rscript 2.fragmentomics_epigenetics/line.R 2.fragmentomics_epigenetics/4.nucleosome/$cohortID.Diff_nuc $CancerType.nuc $CancerType 2.fragmentomics_epigenetics/4.nucleosome/$cohortID.$CancerType.line.pdf
Rscript 2.fragmentomics_epigenetics/line.R 2.fragmentomics_epigenetics/4.nucleosome/$cohortID.Diff_nuc $ControlType.nuc $ControlType 2.fragmentomics_epigenetics/4.nucleosome/$cohortID.$ControlType.line.pdf

Rscript 2.fragmentomics_epigenetics/box.R 2.fragmentomics_epigenetics/4.nucleosome/$cohortID.Diff_nuc Diff_nuc 2.fragmentomics_epigenetics/4.nucleosome/$cohortID.Diff_nuc.box.pdf
Rscript 2.fragmentomics_epigenetics/roc.R 2.fragmentomics_epigenetics/4.nucleosome/$cohortID.Diff_nuc Diff_nuc 2.fragmentomics_epigenetics/4.nucleosome/$cohortID.Diff_nuc.roc.pdf

## E-index
Rscript 2.fragmentomics_epigenetics/line.R 2.fragmentomics_epigenetics/5.Eindex/$cohortID.Diff_Eindex $CancerType.E-index $CancerType 2.fragmentomics_epigenetics/5.Eindex/$cohortID.$CancerType.line.pdf
Rscript 2.fragmentomics_epigenetics/line.R 2.fragmentomics_epigenetics/5.Eindex/$cohortID.Diff_Eindex $ControlType.E-index $ControlType 2.fragmentomics_epigenetics/5.Eindex/$cohortID.$ControlType.line.pdf

Rscript 2.fragmentomics_epigenetics/box.R 2.fragmentomics_epigenetics/5.Eindex/$cohortID.Diff_Eindex Diff_E-index 2.fragmentomics_epigenetics/5.Eindex/$cohortID.Diff_E-index.box.pdf
```

We used the following commands for variant-associated DNA methylation analysis (per sample) in Bie et al. cohort:
```
## Need to update sampleID and path to the FASTQ files
sid=SampleID
FASTQ1=/path/to/$sid.R1.fq.gz
FASTQ2=/path/to/$sid.R2.fq.gz

## extract variants supported by both Wt- and Mut-DNA
perl 2.fragmentomics_epigenetics/6.methylation/extract_snp.pl $sid.sequenced.alleles $sid.common.snv

## extract Wt- and Mut-DNA reads for extracted variants
perl 2.fragmentomics_epigenetics/6.methylation/extract_fragment_id.pl $sid.sequenced.alleles $sid.common.snv $sid
perl 2.fragmentomics_epigenetics/6.methylation/extract_reads.pl $FASTQ1 $FASTQ2 $sid.wt  &
perl 2.fragmentomics_epigenetics/6.methylation/extract_reads.pl $FASTQ1 $FASTQ2 $sid.mut &
wait

## Re-analyze the Wt- and Mut-DNA to get methylation levels
msuite2 -x hg38 -1 $sid.wt.R1.fq.gz -2 $sid.wt.R2.fq.gz -o Msuite2.$sid.wt -k illumina --cut-r1-tail 25 --cut-r2-head 25 --aligner hisat2 -p 16
cd Msuite2.$sid.wt
make && make clean
cd ../

msuite2 -x hg38 -1 $sid.mut.R1.fq.gz -2 $sid.mut.R2.fq.gz -o Msuite2.$sid.mut -k illumina --cut-r1-tail 25 --cut-r2-head 25 --aligner hisat2 -p 16
cd Msuite2.$sid.mut
make && make clean
cd ../
## final methylation calls for Wt- and Mut-DNA were provided in "Processed.files/3.Bie/Wt_DNA_methylation" and "Processed.files/3.Bie/Mut_DNA_methylation" directories

## extract CpG sites in promoters, gene bodies, and intergenic regions
AnnoPath=2.fragmentomics_epigenetics/6.methylation
perl 2.fragmentomics_epigenetics/6.methylation/extract_region.pl Processed.files/3.Bie/Wt_DNA_methylation/$sid.n.CpG.meth.call.gz $AnnoPath $sid.n
perl 2.fragmentomics_epigenetics/6.methylation/extract_region.pl Processed.files/3.Bie/Mut_DNA_methylation/$sid.t.CpG.meth.call.gz $AnnoPath $sid.t

## calculate Diff_methylation in promoters, gene bodies, and intergenic regions
perl 2.fragmentomics_epigenetics/6.methylation/calc_methylation.pl $sid promoter > $sid.promoter.Diff_methylation
perl 2.fragmentomics_epigenetics/6.methylation/calc_methylation.pl $sid genebody > $sid.genebody.Diff_methylation
perl 2.fragmentomics_epigenetics/6.methylation/calc_methylation.pl $sid intergenic > $sid.intergenic.Diff_methylation

## visualize Diff-methylation after analyzing all samples (the Diff-methylation results were provided in this package)
## Here we use genebody as an example.
region=genebody
Rscript 2.fragmentomics_epigenetics/6.methylation/line_methy.R 2.fragmentomics_epigenetics/6.methylation/$region.Diff_methylation $region 2.fragmentomics_epigenetics/6.methylation/$region
Rscript 2.fragmentomics_epigenetics/6.methylation/box_methy.R 2.fragmentomics_epigenetics/6.methylation/$region.Diff_methylation $region 2.fragmentomics_epigenetics/6.methylation/$region
Rscript 2.fragmentomics_epigenetics/6.methylation/roc_methy.R 2.fragmentomics_epigenetics/6.methylation/$region.Diff_methylation $region 2.fragmentomics_epigenetics/6.methylation/$region
```

## 7. Signature deconvolution
For each dataset, we randomly selected 8-10 control samples as "background signatures", and pooled them with COSMIC SBS signatures to deconvolute the rest samples:
```
## COSMIC SBS signatures are downloaded from https://cancer.sanger.ac.uk/signatures/downloads/.
## Note that known sequencing artefacts are NOT used in our analysis.

## The following controls were used as background (recorded in "control_background.txt" files):
## HCC cohort: A67 A46 A60 A74 A96 A77 A97 A89
## Liang cohort: all controls
## Bie cohort: HRR1235773 HRR1235415 HRR1235501 HRR1235571 HRR1235605 HRR1235344 HRR1235464 HRR1235646 HRR1235603 HRR1235527

## The following are commands for analyzing HCC cohort
cohortID=HCC
cut -f 2- 1.variant/$cohortID.control_background.txt | paste 1.variant/COSMIC_signature.txt - > 1.variant/$cohortID.signature.txt
perl 1.variant/select.pl 1.variant/$cohortID.96_mutation_profile.txt 1.variant/$cohortID.control_background.txt 1.variant/$cohortID.96_mutation_profile_final.txt
Rscript 1.variant/deconvolution.R 1.variant/$cohortID.96_mutation_profile_final.txt 1.variant/$cohortID.signature.txt 1.variant/$cohortID

## plot per sample
Rscript 1.variant/Fraction_of_COSMIC_box.R 1.variant/$cohortID.All_SBS_contribution.txt 1.variant/$cohortID.Fraction_of_COSMIC.pdf

## plot per SBS
awk '{for(i=1;i<=NF;i++)a[i]=a[i]?a[i]"\t"$i:$i}END{for(i=1;i<=NF;i++)print a[i]}' 1.variant/$cohortID.COSMIC_SBS_contribution.txt | perl -alne 'print if $F[0] !~ /Ctr/'> 1.variant/$cohortID.tumor_sbs.plot.txt
Rscript 1.variant/SBS_contribution_plot.R 1.variant/$cohortID.tumor_sbs.plot.txt SBS_contribution 1.variant/$cohortID.SBS_contribution.pdf
```

## 8. FreeSV and FreeSV-m models
FreeSV model was built through integrating genomic, fragmentomic, and epigenetic features associated with somatic variants in Bie et al. cohort. We pooled the individual features calculated by above commands, and sorted them in "3.AI_model/FreeSV_data.txt" file. We used the following commands to build and evaluate the model:
```
## build model
Rscript 3.AI_model/GBM_parallel.R 3.AI_model/FreeSV_data.txt FreeSV_test
perl 3.AI_model/mean_pred.pl FreeSV_test > FreeSV_test.pred.txt

## ROC
Rscript 3.AI_model/roc.R 3.AI_model/FreeSV_test.pred.txt Test 3.AI_model/FreeSV_test.roc.pdf
Rscript 3.AI_model/roc_stage.R 3.AI_model/FreeSV_test.pred.txt Test_stage 3.AI_model/FreeSV_test_stage.roc.pdf
Rscript 3.AI_model/box.R 3.AI_model/FreeSV_test.pred.txt Test_pred 3.AI_model/FreeSV_pred.pdf
Rscript 3.AI_model/box_stage.R 3.AI_model/FreeSV_test.pred.txt Test_pred_stage 3.AI_model/FreeSV_pred_stage.pdf
```

For FreeSV-m model, we pooled the features used in FreeSV and 4 genomewide features calculated in Bie et al. work. We then built the model using the "training" group and applied the model on "testing" group (both defined in Bie et al. study) to evaluate its performance. The feature matrices were recorded in "FreeSV-m_train_data.txt" and "3.AI_model/FreeSV-m_test_data.txt" files. We used the following commands to build and evaluate the model:
```
## train model on training group
Rscript 3.AI_model/GBM_parallel.R 3.AI_model/FreeSV-m_train_data.txt FreeSV-m_train
perl 3.AI_model/mean_pred.pl FreeSV-m_train > FreeSV-m_train.pred.txt

## apply FreeSV-m model on testing group
Rscript 3.AI_model/apply_model.R 3.AI_model/FreeSV-m_test_data.txt FreeSV-m_train 3.AI_model/FreeSV-m_test
perl 3.AI_model/mean_pred.pl FreeSV-m_test > FreeSV-m_test.pred.txt

## ROC on training group
Rscript 3.AI_model/roc.R 3.AI_model/FreeSV-m_train.pred.txt Train 3.AI_model/FreeSV-m_train.roc.pdf
Rscript 3.AI_model/roc_stage.R 3.AI_model/FreeSV-m_train.pred.txt Train_stage 3.AI_model/FreeSV-m_train_stage.roc.pdf
Rscript 3.AI_model/box.R 3.AI_model/FreeSV-m_train.pred.txt Train_pred 3.AI_model/FreeSV-m_train_pred.pdf
Rscript 3.AI_model/box_stage.R 3.AI_model/FreeSV-m_train.pred.txt Train_pred_stage 3.AI_model/FreeSV-m_train_pred_stage.pdf

## ROC on testing group
Rscript 3.AI_model/roc.R 3.AI_model/FreeSV-m_test.pred.txt Test 3.AI_model/FreeSV-m_test.roc.pdf
Rscript 3.AI_model/roc_stage.R 3.AI_model/FreeSV-m_test.pred.txt Test_stage 3.AI_model/FreeSV-m_test_stage.roc.pdf
Rscript 3.AI_model/box.R 3.AI_model/FreeSV-m_test.pred.txt Test_pred 3.AI_model/FreeSV-m_test_pred.pdf
Rscript 3.AI_model/box_stage.R 3.AI_model/FreeSV-m_test.pred.txt Test_pred_stage 3.AI_model/FreeSV-m_test_pred_stage.pdf
```

## 9. Murine cfDNA data
We identified the somatic variants, and analyzed the cfDNA size and end motif usages associated with Wt- and Mut-DNA using the following commands:
```
## Need to update sampleID and path to the FASTQ files
sid=SampleID
FASTQ1=/path/to/$sid.R1.fq.gz
FASTQ2=/path/to/$sid.R2.fq.gz

## data preprocess
ktrim -1 $FASTQ1 -2 $FASTQ2 -t 8 -o $sid.ktrim -k BGI

## read alignment. Need to build C57BL/6J genome index for BWA first and replace the following path
## the genome fasta file is downloaded from ftp://ftp.jax.org/b6eve
C57index=/path/to/C57.bwa.index
bwa mem -t 16 -M -Y $C57index $sid.ktrim.read1.fq $sid.ktrim.read2.fq | samtools view -b -@ 16 - | samtools sort -@ 16 -o $sid.sort.bam - 2> $sid.bwa.log

## remove duplicates
java -Xmx16g -XX:ParallelGCThreads=8 -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$sid.sort.bam O=$sid.mkdup.bam METRICS_FILE=$sid.mkdup.bam.mat TMP_DIR=tmp 2>$sid.mkdup.log

## call variants. Need to update the path to C57BL/6J genome file (FASTA format)
C57fasta=/path/to/C57.fa
samtools view -@ 4 -h $sid.mkdup.bam | perl -alne 'print and next if /^@/; next if $F[1] & 0x100 || $F[2] eq "*"; print if ($F[1] & 0x02) && ($F[6] eq "=")' | samtools view -b -@ 16 -o $sid.filter.bam
samtools index -@ 16 $sid.filter.bam
bcftools mpileup -q 60 -Q 30 -Ov -f $C57fasta --threads 16 $sid.filter.bam | bcftools call -vc -V indels -Oz --threads 16 -o $sid.vcf.gz

## filter variants
zcat $sid.vcf.gz | perl -alne 'print and next if /^#/; print if $F[0]=~/^chr\d+$/ && $F[5] >= 30 && $F[4] !~ /,/' | gzip > $sid.hdfilt.pass.snp.raw.vcf.gz
zcat $sid.hdfilt.pass.snp.raw.vcf.gz | perl -alne 'next if /^#/; print "$F[0]\t$F[1]"' | perl 1.variant/get_dmbs_pos.pl - | uniq > $sid.dmbs.txt
perl 1.variant/filt_vcf_dmbs.pl $sid.hdfilt.pass.snp.raw.vcf.gz $sid.dmbs.txt $sid.hdfilt.pass.snp.tmp.vcf.gz
perl 1.variant/sub_dp_mouse.pl  $sid.hdfilt.pass.snp.tmp.vcf.gz $sid.final.vcf.gz

## extract Wt- and Mut-DNA; the processed files are under Processed.files/4.Mouse/ directory.
WtDNA=Processed.files/4.Mouse/Wt_DNA/$sid.bed.gz
MutDNA=Processed.files/4.Mouse/Mut_DNA/$sid.bed.gz

## size
perl 2.fragmentomics_epigenetics/2.size/bed2size.pl $WtDNA  $sid.Wt.size
perl 2.fragmentomics_epigenetics/2.size/bed2size.pl $MutDNA $sid.Mut.size

Rscript 2.fragmentomics_epigenetics/2.size/plot_size.R $sid
n=`perl 2.fragmentomics_epigenetics/2.size/size_freq.pl $sid.Wt.size`
t=`perl 2.fragmentomics_epigenetics/2.size/size_freq.pl $sid.Mut.size`
echo -e "$sid\t$n\t$t" | perl -alne 'print join "\t", $F[0], $F[1], $F[2], $F[2]-$F[1]' > $sid.Diff_size

## end motif
perl 2.fragmentomics_epigenetics/3.motif/grab.end.with.extension.pl $WtDNA  $sid.Wt  PE &
perl 2.fragmentomics_epigenetics/3.motif/grab.end.with.extension.pl $MutDNA $sid.Mut PE &
wait
perl 2.fragmentomics_epigenetics/3.motif/bed2fa.local.pl $C57fasta $sid.Wt.left.outer2.inner4.bed $sid.Mut.left.outer2.inner4.bed
perl 2.fragmentomics_epigenetics/3.motif/extract.motif.with.extension.pl $sid.Wt.left.outer2.inner4.fa  > $sid.Wt.motif  &
perl 2.fragmentomics_epigenetics/3.motif/extract.motif.with.extension.pl $sid.Mut.left.outer2.inner4.fa > $sid.Mut.motif &
wait

ccca_n=$(grep -w "^CCCA" $sid.Wt.motif  | cut -f 5)
ccca_t=$(grep -w "^CCCA" $sid.Mut.motif | cut -f 5)
echo -e "$sid\t$ccca_n\t$ccca_t" | perl -alne 'print join "\t", $F[0], $F[1], $F[2], $F[2]-$F[1]' > $sid.motif.ccca

ctcc_n=$(grep -w "^CTCC" $sid.Wt.motif  | cut -f 3)
ctcc_t=$(grep -w "^CTCC" $sid.Mut.motif | cut -f 3)
echo -e "$sid\t$ctcc_n\t$ctcc_t" | perl -alne 'print join "\t", $F[0], $F[1], $F[2], $F[2]-$F[1]' > $sid.motif.ctcc
```
