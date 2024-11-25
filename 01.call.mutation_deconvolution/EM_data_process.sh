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
  
  /usr/bin/msuite2 -x hg38 -1 data/$sid.1.fq.gz -2 data/$sid.2.fq.gz -o $sid -k illumina --cut-r1-tail 25 --cut-r2-head 25 --aligner hisat2

  cd $sid

  make && make clean

  cd ../

  echo "Sample $sid processing complete."
done < "$sample_list"

