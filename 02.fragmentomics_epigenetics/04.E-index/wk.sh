#!/bin/bash
#
# For calculating E-index values in cfDNA, please refer to An et al. Nature Communications 2023; 10.1038/s41467-023-35959-6. [PubMed link] (https://pubmed.ncbi.nlm.nih.gov/36653380/).
#

control_e="control.e"
hcc_e="hcc.e"
diff_txt="diff.txt"

echo -ne "Sample\tWt\tMut\n" > $control_e
echo -ne "Sample\tWt\tMut\n" > $hcc_e

paste n.txt t.txt | cut -f 1,5,10 | head -n 24 >> $control_e
paste n.txt t.txt | cut -f 1,5,10 | tail -n 56 >> $hcc_e

Rscript line.R \
  $control_e Control_E_index control.pdf

Rscript line.R \
  $hcc_e HCC_E_index hcc.pdf

# Diff-E-index
echo -ne "Sample\tGroup\tDiff\n" > $diff_txt

sed 1d $hcc_e | perl -alne 'print join "\t", $F[0], "Tumor", $F[2]-$F[1]' >> $diff_txt
sed 1d $control_e | perl -alne 'print join "\t", $F[0], "Control", $F[2]-$F[1]' >> $diff_txt

Rscript box.R \
  $diff_txt E_index_diff diff.pdf

