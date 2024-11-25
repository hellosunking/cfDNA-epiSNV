#!/bin/bash

sample_hcc="sample.hcc"
sample_control="sample.control"

output_dir="../02.size"
mkdir -p $output_dir

hcc_output="${output_dir}/hcc.size"
echo -ne "Sample\tWt\tMut\n" > $hcc_output
for i in `cut -f 1 $sample_hcc`; do
  a=$(perl size_freq.pl "$output_dir/$i.tumor.size")
  b=$(perl size_freq.pl "$output_dir/$i.normal.size")
  echo -ne "$i\t$b\t$a\n" >> $hcc_output
done

control_output="${output_dir}/control.size"
echo -ne "Sample\tWt\tMut\n" > $control_output
for i in `cut -f 1 $sample_control`; do
  a=$(perl size_freq.pl "$output_dir/$i.tumor.size")
  b=$(perl size_freq.pl "$output_dir/$i.normal.size")
  echo -ne "$i\t$b\t$a\n" >> $control_output
done

Rscript_path="/lustre/home/zyzhang/04.software/biosoft/miniconda/miniconda3/bin/Rscript"
Rscript line.R $hcc_output "HCC_S150_cumulation" "${output_dir}/hcc.pdf"
Rscript line.R $control_output "Control_S150_cumulation" "${output_dir}/control.pdf"

# Diff-Size
diff_output="${output_dir}/diff.txt"
echo -ne "Sample\tGroup\tDiff\n" > $diff_output
sed 1d $hcc_output | perl -alne 'print join "\t", $F[0], "Tumor", $F[2]-$F[1]' >> $diff_output
sed 1d $control_output | perl -alne 'print join "\t", $F[0], "Control", $F[2]-$F[1]' >> $diff_output

Rscript box.R $diff_output "S150_percent_diff" "${output_dir}/diff.pdf"

Rscript roc.pdf
