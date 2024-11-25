library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

args <- commandArgs(TRUE)
input <- args[1]
prefix     <- args[2]

dat <- read.table(file = input, header = TRUE, check.names = FALSE, row.names = 1)
dat <- t(dat)
dat <- as.matrix(dat)
dat <- scale(dat)

dat_transposed <- t(dat)

color_pal <- brewer.pal(n = 7, name = "RdYlBu")
color_pal <- rev(color_pal)

column_colors <- rep("black", ncol(dat_transposed))
column_colors[grepl("HCC", colnames(dat_transposed))] <- "#CD2626"
column_colors[grepl("Ctr", colnames(dat_transposed))] <- "grey"

set.seed(7)
km_res <- kmeans(dat, centers = 2)
cluster <- km_res$cluster

column_anno <- HeatmapAnnotation(
  Cluster = factor(cluster),
  col = list(Cluster = c("1" = "grey", "2" = "#CD2626"))
)

pdf(file= paste(prefix, ".pdf", sep="") )
Heatmap(dat_transposed, name = "Profile",
        col = color_pal,
        column_names_gp = gpar(col = column_colors, fontsize = 4),
        row_names_gp = gpar(fontsize = 4),
        show_column_names = TRUE,
        show_row_names = T,
        set.seed(123),
        column_km = 2,
        width = unit(15, "cm"),
        height = unit(12, "cm"),
        column_dend_side = "top",
        column_dend_height = unit(1, "cm"),
        top_annotation = column_anno
)
dev.off()
