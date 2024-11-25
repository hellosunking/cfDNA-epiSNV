#!/bin/bash
#
# For calculating nucleosome track in cfDNA, please refer to An et al. Nature Communications 2023; 10.1038/s41467-023-35959-6. [PubMed link] (https://pubmed.ncbi.nlm.nih.gov/36653380/).
#

output_dir="analysis"
mkdir -p $output_dir

for t in `cut -f 1 sample.control`
do
  a=$(perl -alne 'print $F[1] if $F[0] >= -50 && $F[0] <= 50' ctr_mut/$t.all.D.dist | awk '{sum+=$1}END{print sum}')
  b=$(awk '{sum+=$2}END{print sum}' ctr_mut/$t.all.D.dist)
  c=$(perl -alne 'print $F[1] if $F[0] >= -50 && $F[0] <= 50' ctr_mut/$t.all.U.dist | awk '{sum+=$1}END{print sum}')
  d=$(awk '{sum+=$2}END{print sum}' ctr_mut/$t.all.U.dist)

  a1=$(perl -alne 'print $F[1] if $F[0] >= -50 && $F[0] <= 50' ctr_wt/$t.all.D.dist | awk '{sum+=$1}END{print sum}')
  b1=$(awk '{sum+=$2}END{print sum}' ctr_wt/$t.all.D.dist)
  c1=$(perl -alne 'print $F[1] if $F[0] >= -50 && $F[0] <= 50' ctr_wt/$t.all.U.dist | awk '{sum+=$1}END{print sum}')
  d1=$(awk '{sum+=$2}END{print sum}' ctr_wt/$t.all.U.dist)

  echo -ne "$t\t$a\t$b\t$a1\t$b1\n" >> $output_dir/ctr.d.tmp50
  echo -ne "$t\t$c\t$d\t$c1\t$d1\n" >> $output_dir/ctr.u.tmp50
done

for t in `cut -f 1 sample.hcc`
do
  a=$(perl -alne 'print $F[1] if $F[0] >= -50 && $F[0] <= 50' hcc_mut/$t.all.D.dist | awk '{sum+=$1}END{print sum}')
  b=$(awk '{sum+=$2}END{print sum}' hcc_mut/$t.all.D.dist)
  c=$(perl -alne 'print $F[1] if $F[0] >= -50 && $F[0] <= 50' hcc_mut/$t.all.U.dist | awk '{sum+=$1}END{print sum}')
  d=$(awk '{sum+=$2}END{print sum}' hcc_mut/$t.all.U.dist)

  a1=$(perl -alne 'print $F[1] if $F[0] >= -50 && $F[0] <= 50' hcc_wt/$t.all.D.dist | awk '{sum+=$1}END{print sum}')
  b1=$(awk '{sum+=$2}END{print sum}' hcc_wt/$t.all.D.dist)
  c1=$(perl -alne 'print $F[1] if $F[0] >= -50 && $F[0] <= 50' hcc_wt/$t.all.U.dist | awk '{sum+=$1}END{print sum}')
  d1=$(awk '{sum+=$2}END{print sum}' hcc_wt/$t.all.U.dist)

  echo -ne "$t\t$a\t$b\t$a1\t$b1\n" >> $output_dir/hcc.d.tmp50
  echo -ne "$t\t$c\t$d\t$c1\t$d1\n" >> $output_dir/hcc.u.tmp50
done

paste $output_dir/ctr.d.tmp50 $output_dir/ctr.u.tmp50 | \
perl -alne 'print join "\t", $F[0], $F[1]+$F[6], $F[2]+$F[7], $F[3]+$F[8], $F[4]+$F[9]' > $output_dir/ctr.all.tmp50

paste $output_dir/hcc.d.tmp50 $output_dir/hcc.u.tmp50 | \
perl -alne 'print join "\t", $F[0], $F[1]+$F[6], $F[2]+$F[7], $F[3]+$F[8], $F[4]+$F[9]' > $output_dir/hcc.all.tmp50

perl -alne 'print join "\t", "Sample", "Mut", "Wt" if $.==1; print join "\t", $F[0], $F[1]/$F[2], $F[3]/$F[4]' \
  $output_dir/ctr.all.tmp50 > $output_dir/ctr.all.txt50

perl -alne 'print join "\t", "Sample", "Mut", "Wt" if $.==1; print join "\t", $F[0], $F[1]/$F[2], $F[3]/$F[4]' \
  $output_dir/hcc.all.tmp50 > $output_dir/hcc.all.txt50

Rscript line.R \
  $output_dir/hcc.all.txt50 HCC_U_D_End $output_dir/hcc.all50_new.pdf

Rscript line.R \
  $output_dir/ctr.all.txt50 Control_U_D_End $output_dir/ctr.all50_new.pdf

# Diff-nucleosome
echo -ne "Sample\tGroup\tDiff\n" > $output_dir/all.diff.txt50

sed 1d $output_dir/hcc.all.txt50 | \
perl -alne 'print join "\t", $F[0], "Tumor", $F[1]-$F[2]' >> $output_dir/all.diff.txt50

sed 1d $output_dir/ctr.all.txt50 | \
perl -alne 'print join "\t", $F[0], "Control", $F[1]-$F[2]' >> $output_dir/all.diff.txt50

Rscript box.R \
  $output_dir/all.diff.txt50 U_D_End_percent_diff $output_dir/all.diff50_new.pdf

Rscript roc.R
