library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript script.R <input> <title> <output_prefix>\n")
}
 
input_file <- args[1]
plot_title <- args[2]
output <- args[3]

data <- read.delim(input_file, sep = "\t")

output_pdf <- paste0(output, ".methylation.box.pdf")

title_final <- paste0("Methylation.", plot_title)

data$Group <- factor(data$Group, levels = c("Control", "BRCA", "COREAD", "ESCA", "STAD", "LIHC", "NSCLC", "PACA"))

colors <- c("Control" = "grey50", "LIHC" = "#CD2626", "BRCA" = "darkorange", "COREAD" = "blue", "ESCA" = "deepskyblue", "STAD" = "turquoise4", "NSCLC" = "purple4", "PACA" = "darkorange4")

p <- ggplot(data, aes(x = Group, y = Diff*100, color = Group)) +
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.6, color = "grey") +
  geom_jitter(width = 0.2, size = 0.8) +
  scale_color_manual(values = colors) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
	axis.ticks.x = element_line(color = "black", size = 0.2),
	axis.ticks.y = element_line(color = "black", size = 0.2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5)) +
  labs(x = "Group", y = "Diff_percent(%)", title = title_final) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, comparisons = list(c("Control", "BRCA"), c("Control", "COREAD"), c("Control", "ESCA"), c("Control", "STAD"), c("Control", "LIHC"), c("Control", "NSCLC"), c("Control", "PACA")), label = "p.format", p.format = function(x) format(x, digitis = 1))

ggsave(output_pdf, plot = p, width = 6, height = 6)
