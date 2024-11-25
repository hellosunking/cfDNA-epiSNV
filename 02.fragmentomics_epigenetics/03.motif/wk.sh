#!/bin/bash
#
# For calculating end motif in cfDNA, please refer to Jiang et al. Cancer Discovery 2020 May; 10(5):664-673. [PubMed link] (https://www.ncbi.nlm.nih.gov/pubmed/32111602).
#

sample_hcc="sample.hcc"
sample_control="sample.control"

motif_dir="../03.motif"
mkdir -p $motif_dir

echo -ne "Sample\tWt\tMut\n" > $motif_dir/hcc.motif.ccca
echo -ne "Sample\tWt\tMut\n" > $motif_dir/hcc.motif.ctcc

for i in $(cut -f 1 $sample_hcc); do
  ccca_tumor=$(grep -w "CCCA" $motif_dir/$i.tumor.left.motif | cut -f 5)
  ccca_normal=$(grep -w "CCCA" $motif_dir/$i.normal.left.motif | cut -f 5)
  ctcc_tumor=$(grep -w "CTCC" $motif_dir/$i.tumor.left.motif | cut -f 3)
  ctcc_normal=$(grep -w "CTCC" $motif_dir/$i.normal.left.motif | cut -f 3)

  echo -ne "$i\t$ccca_normal\t$ccca_tumor\n" >> $motif_dir/hcc.motif.ccca
  echo -ne "$i\t$ctcc_normal\t$ctcc_tumor\n" >> $motif_dir/hcc.motif.ctcc
done

echo -ne "Sample\tWt\tMut\n" > $motif_dir/control.motif.ccca
echo -ne "Sample\tWt\tMut\n" > $motif_dir/control.motif.ctcc

for i in $(cut -f 1 $sample_control); do
  ccca_tumor=$(grep -w "CCCA" $motif_dir/$i.tumor.left.motif | cut -f 5)
  ccca_normal=$(grep -w "CCCA" $motif_dir/$i.normal.left.motif | cut -f 5)
  ctcc_tumor=$(grep -w "CTCC" $motif_dir/$i.tumor.left.motif | cut -f 3)
  ctcc_normal=$(grep -w "CTCC" $motif_dir/$i.normal.left.motif | cut -f 3)

  echo -ne "$i\t$ccca_normal\t$ccca_tumor\n" >> $motif_dir/control.motif.ccca
  echo -ne "$i\t$ctcc_normal\t$ctcc_tumor\n" >> $motif_dir/control.motif.ctcc
done

Rscript line.R \
  $motif_dir/hcc.motif.ccca HCC_CCCA $motif_dir/hcc.motif.ccca.pdf

Rscript line.R \
  $motif_dir/hcc.motif.ctcc HCC_CTCC $motif_dir/hcc.motif.ctcc.pdf

Rscript line.R \
  $motif_dir/control.motif.ccca Control_CCCA $motif_dir/control.motif.ccca.pdf

Rscript line.R \
  $motif_dir/control.motif.ctcc Control_CTCC $motif_dir/control.motif.ctcc.pdf

echo -ne "Sample\tGroup\tDiff\n" > $motif_dir/ccca.diff.txt
sed 1d $motif_dir/hcc.motif.ccca | perl -alne 'print join "\t", $F[0], "Tumor", $F[2]-$F[1]' >> $motif_dir/ccca.diff.txt
sed 1d $motif_dir/control.motif.ccca | perl -alne 'print join "\t", $F[0], "Control", $F[2]-$F[1]' >> $motif_dir/ccca.diff.txt

echo -ne "Sample\tGroup\tDiff\n" > $motif_dir/ctcc.diff.txt
sed 1d $motif_dir/hcc.motif.ctcc | perl -alne 'print join "\t", $F[0], "Tumor", $F[2]-$F[1]' >> $motif_dir/ctcc.diff.txt
sed 1d $motif_dir/control.motif.ctcc | perl -alne 'print join "\t", $F[0], "Control", $F[2]-$F[1]' >> $motif_dir/ctcc.diff.txt

# Diff-Motif
Rscript box_motif.R \
  $motif_dir/ccca.diff.txt CCCA_percent_diff $motif_dir/ccca.diff.pdf

Rscript box_motif.R \
  $motif_dir/ctcc.diff.txt CTCC_percent_diff $motif_dir/ctcc.diff.pdf

Rscript roc.R
