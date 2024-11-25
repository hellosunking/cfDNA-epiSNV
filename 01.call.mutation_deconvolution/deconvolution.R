library(MutationalPatterns)
library("gridExtra")
library(BSgenome)
library("NMF")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

args <- commandArgs(TRUE)
input <- args[1] 
signature <- args[2]

mut_mat <- read.table(input, header = TRUE, row.names = 1)
signatures = read.table(signature, sep = "\t", header = TRUE, row.names = 1)

selected_columns <- sample(1:24, 8)
selected_mut_mat <- mut_mat[, selected_columns, drop = FALSE]

column_sums <- colSums(selected_mut_mat)

selected_mut_mat <- selected_mut_mat / column_sums

signatures <- cbind(signatures, selected_mut_mat)

mut_mat <- mut_mat[, -selected_columns, drop = FALSE]

meta_cols <- c(1)
signatures <- as.matrix(signatures[, -meta_cols, drop = FALSE])
fit_res <- fit_to_signatures(mut_mat, signatures)

filtered_heatmap_matrix <- fit_res$contribution

write.table(filtered_heatmap_matrix, paste("all_refsig_sample_contribution_count.txt",sep=""), sep = "\t", quote = FALSE)
filtered_heatmap_matrix <- t(filtered_heatmap_matrix)
# relative contribution
filtered_heatmap_matrix_norm <- filtered_heatmap_matrix / rowSums(filtered_heatmap_matrix)
filtered_heatmap_matrix_norm <- t(filtered_heatmap_matrix_norm)
filtered_heatmap_matrix <- t(filtered_heatmap_matrix)
write.table(filtered_heatmap_matrix_norm, paste("all_refsig_sample_contribution_rate.txt",sep=""), sep = "\t", quote = FALSE)

new_variable <- filtered_heatmap_matrix_norm[!grepl("Ctr", rownames(filtered_heatmap_matrix_norm)), ]
write.table(new_variable, paste("SBS_sample_contribution_rate.txt",sep=""), sep = "\t", quote = FALSE)

pdf("all_ref_sig_oriprof_cossi.pdf", width = 12, height = 15)
plot_original_vs_reconstructed(mut_mat, fit_res$reconstructed, y_intercept = 0.95)
dev.off()

