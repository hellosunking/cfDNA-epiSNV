#!/bin/bash
#
# This program is designed for AI-model FreeSV
#

data_file=$1
sample_info_file=$2

# Check if the input files exist
if [[ ! -f $data_file ]]; then
  echo "Error: Data file '$data_file' not found!"
  exit 1
fi

if [[ ! -f $sample_info_file ]]; then
  echo "Error: Sample info file '$sample_info_file' not found!"
  exit 1
fi

# Step 1: Run the R script with data file
Rscript GBM_parallel.R "$data_file"

# Step 2: Generate predictions
for i in `seq 1 100`; do
  for j in `seq 1 10`; do
    sed 1d pred/test_pred_rep_"$i"_fold_"$j".txt | cut -f 1,3 >> test.tmp
    sed 1d pred/train_pred_rep_"$i"_fold_"$j".txt | cut -f 1,3 >> train.tmp
  done
done

# Step 3: Calculate mean predictions
perl mean_pred.pl test.tmp test.tmp1
perl mean_pred.pl train.tmp train.tmp1

# Step 4: Clean temporary files
rm test.tmp train.tmp

# Step 5: Generate headers and row data
echo -ne "Sid\tGroup\tpred\n" > head
cut -f 1,2 "$sample_info_file" > row

# Step 6: Combine predictions with sample info
paste row test.tmp1 | cut -f 1,2,4 > test.tmp2
paste row train.tmp1 | cut -f 1,2,4 > train.tmp2

# Step 7: Clean intermediate files
rm test.tmp1 train.tmp1

# Step 8: Create final prediction files
cat head test.tmp2 > test.pred.txt
cat head train.tmp2 > train.pred.txt

rm test.tmp2 train.tmp2 head row

# Step 9: Generate ROC plot
Rscript roc.R test.pred.txt Train train.roc.pdf
Rscript roc.R test.pred.txt Test test.roc.pdf

# Step 10: Prediction scores
Rscript pred_scores_box.R pred.scores.txt Pred pred.pdf
