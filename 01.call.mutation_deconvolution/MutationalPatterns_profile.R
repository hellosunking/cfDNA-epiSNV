input_file <- commandArgs(trailingOnly = TRUE)[1]
prefix <- commandArgs(trailingOnly = TRUE)[2]

library(MutationalPatterns)
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

vcfinfo <- read.table(input_file, header = FALSE, stringsAsFactors = FALSE)
colnames(vcfinfo) <- c("sample", "group", "vcf_path")

grl <- read_vcfs_as_granges(vcfinfo$vcf_path, vcfinfo$sample, ref_genome, type = "all")
tissue <- vcfinfo$group

snv_grl <- get_mut_type(grl, type = "snv")
type_occurrences <- mut_type_occurrences(snv_grl, ref_genome)

pdf(file= paste(prefix, ".sper_plots.pdf", sep=""), width = 10, height = 10)
p1 <- plot_spectrum(type_occurrences)
p4 <- plot_spectrum(type_occurrences, by = tissue, CT = TRUE, legend = TRUE)
print(p1)
print(p4)
dev.off()

library("gridExtra")
pdf(file= paste(prefix, ".output_plots.pdf", sep=""), width = 15, height = 10)
grid.arrange(p1, p4, ncol = 2, widths = c(1, 1))
dev.off()

mut_mat <- mut_matrix(vcf_list=snv_grl, ref_genome = ref_genome)
pdf(file= paste(prefix, ".96_profile_plot.pdf", sep=""), width = 15, height = 150)
plot_96_profile(mut_mat)
dev.off()

write.table(type_occurrences, paste(prefix,".type_occurrences.xls",sep=""), sep = "\t", quote = FALSE)
write.table(mut_mat, paste(prefix,".mut_mat.xls",sep=""), sep = "\t", quote = FALSE)
