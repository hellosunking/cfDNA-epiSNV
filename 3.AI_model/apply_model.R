library(caret)
library(pROC)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript script.R <input> <model_prefix> <output_prefix>\n")
}

input_file <- args[1]      
prefix <- args[2]          
output <- args[3]     

test_data <- read.table(input_file, header = TRUE, row.names = 1)

test_Type <- test_data$Type

for (rep in 1:100) {
  for (fold in 1:10) {
    
    model_file <- paste0(prefix, ".rep_", rep, "_fold_", fold, ".rds")
    
    model <- readRDS(model_file)
    
    test_pred <- predict(model, newdata = test_data, type = "prob")
    
    test_pred_t <- data.frame(sid = rownames(test_data),
                              Type = test_Type,
                              pred_class1 = test_pred$class1)
    
    output_file_name <- paste0(output, ".test.rep_", rep, "_fold_", fold, ".txt")
    
    write.table(test_pred_t, file = output_file_name, sep = "\t", col.names = F, row.names = FALSE, quote = FALSE)
  }
}
