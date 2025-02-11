args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript script.R <input_file> <output_file>\n")
}

input_file <- args[1]
output <- args[2]

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

mut_mat <- mut_matrix(vcf_list=snv_grl, ref_genome = ref_genome)

write.table(mut_mat, output, sep = "\t", quote = FALSE)
