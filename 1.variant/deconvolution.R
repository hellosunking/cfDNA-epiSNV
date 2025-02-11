library(MutationalPatterns)
library("gridExtra")
library(BSgenome)
library("NMF")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

args <- commandArgs(TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript script.R <profile> <signature> <output_prefix>\n")
}

input <- args[1]
sig <- args[2]
prefix <- args[3]

mut_mat <- read.table(input, header = TRUE, row.names = 1)
rownames(mut_mat)[1:16] <- paste("Ctr_", rownames(mut_mat)[1:16], sep = "")
rownames(mut_mat)[17:nrow(mut_mat)] <- paste("HCC_", rownames(mut_mat)[17:nrow(mut_mat)], sep = "")
mut_mat <- t(mut_mat)
signatures = read.table(sig, sep = "\t", header = TRUE)
meta_cols <- c(1)
signatures <- as.matrix(signatures[, -meta_cols, drop = FALSE])
fit_res <- fit_to_signatures(mut_mat, signatures)

filtered_heatmap_matrix <- fit_res$contribution

filtered_heatmap_matrix <- t(filtered_heatmap_matrix)
# relative contribution
filtered_heatmap_matrix_norm <- filtered_heatmap_matrix / rowSums(filtered_heatmap_matrix)
filtered_heatmap_matrix_norm <- t(filtered_heatmap_matrix_norm)
filtered_heatmap_matrix <- t(filtered_heatmap_matrix)
write.table(filtered_heatmap_matrix_norm, paste(prefix, ".All_SBS_contribution.txt",sep=""), sep = "\t", quote = FALSE)

new_variable <- filtered_heatmap_matrix_norm[!grepl("A", rownames(filtered_heatmap_matrix_norm)), ]
write.table(new_variable, paste(prefix, ".COSMIC_SBS_contribution.txt",sep=""), sep = "\t", quote = FALSE)
