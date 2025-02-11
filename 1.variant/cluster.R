library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

args <- commandArgs(TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript script.R <input_file> <sample_info> <output_prefix>\n")
}

input.file <- args[1]
class.file <- args[2]
prefix     <- args[3]

dat <- read.table(file = input.file, header = TRUE, check.names = FALSE, row.names = 1)
dat <- as.matrix(dat)
dat <- scale(dat)

sample_info <- read.table(class.file, header = FALSE, stringsAsFactors = FALSE)
colnames(sample_info) <- c("sample", "group", "other")
sample_info <- sample_info[sample_info$group %in% c("Control", "HCC"), ]

dat <- dat[, rownames(dat) %in% sample_info$sample]

# Transpose data matrix for clustering
dat_transposed <- t(dat)

# Define color palette for heatmap
color_pal <- brewer.pal(n = 7, name = "RdYlBu")
color_pal <- rev(color_pal)

# Initialize column_colors
column_colors <- rep("black", ncol(dat_transposed))  # Initialize with black color

# Assign colors based on group
column_colors[sample_info$group == "Control"] <- "grey"
column_colors[sample_info$group == "HCC"] <- "#CD2626"

# Perform k-means clustering
set.seed(7)
km_res <- kmeans(dat, centers = 2)
cluster <- km_res$cluster

# Create HeatmapAnnotation with group information from sample_info
column_anno <- HeatmapAnnotation(
  Cluster = factor(cluster),
  SampleType = anno_simple(sample_info$group, col = c("Control" = "grey", "HCC" = "#CD2626")),
  col = list(Cluster = c("1" = "#dae3f3", "2" = "#fbe5d6"))
)

# Custom legend for SampleType with specified colors
sampletype_legend <- Legend(
  labels = c("Control", "HCC"), 
  legend_gp = gpar(fill = c("grey", "#CD2626")), 
  title = "Sample Type"
)

# Create heatmap with selected features highlighted
heatmap <- Heatmap(dat_transposed, name = "Profile",
                   col = color_pal,
                   column_names_gp = gpar(col = column_colors, fontsize = 4),
                   row_names_gp = gpar(fontsize = 5.7),
                   show_column_names = FALSE,  # Do not show column names
                   show_row_names = F,
                   set.seed(7),
                   cluster_columns = TRUE,
                   column_km = 2,
                   width = unit(18, "cm"),  # Increase width for wider cells
                   height = unit(18, "cm"),  # Adjust height as needed
                   column_dend_side = "top",
                   column_dend_height = unit(1, "cm"),
                   top_annotation = column_anno
)

draw(heatmap, annotation_legend_list = list(sampletype_legend))

output_file <- paste0(prefix, ".heatmap.pdf")

ggsave(output_file, height = 8, width = 8)
