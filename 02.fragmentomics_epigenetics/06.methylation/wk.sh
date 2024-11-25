#!/bin/bash
# This program is designed for EM-seq data analysis only
#

sample=$1

input_dir="./"
output_dir="./common"
script_dir="./"
mkdir -p $output_dir

# Step1: get variant loci with both Wt-DNA and Mut-DNA
cut -f 1,4,7 $input_dir/$sample.sequenced.alleles | \
sed 's/,/\t/g' | \
perl -alne 'print join "_", $F[0], @F[2..4] if $F[1] eq $F[3]' | \
sort | uniq > $input_dir/$sample.raw.n

cut -f 1,4,7 $input_dir/$sample.sequenced.alleles | \
sed 's/,/\t/g' | \
perl -alne 'print join "_", $F[0], @F[2..4] if $F[1] eq $F[4]' | \
sort | uniq > $input_dir/$sample.raw.t

cat $input_dir/$sample.raw.n $input_dir/$sample.raw.t | \
sort | uniq -d > $input_dir/$sample.snp.id

perl sub_new.pl $input_dir/$sample.sequenced.alleles $input_dir/$sample.snp.id | \
sort -k1,1V -k2,2n > $output_dir/$sample.sequenced.alleles

# Step2: extract the Wt- and Mut-DNA && calculate methylation levles using Msuite2 respectively
# Manually generate files containing the methylation levels of Wt- and Mut-DNA respectively like fragmentomic features analysis  
# Sample	Wt	Mut
Rscript line.R wt_mut.txt Methylation.pdf

# Step3: Diff-Methylation
echo -ne "Sample\tGroup\tDiff\n" > methy.diff.txt
sed 1d wt_mut.control.txt | perl -alne 'print join "\t", $F[0], "Control", $F[2]-$F[1]' >> methy.diff.txt
sed 1d wt_mut.tumor.txt | perl -alne 'print join "\t", $F[0], "Tumor", $F[2]-$F[1]' >> methy.diff.txt

Rscript box.R methy.diff.txt Diff-Methylation diff.pdf

Rscript roc.R methy.diff.txt Diff-Methylation roc.pdf
