library(ggplot2)
library(dplyr)
library(ggpubr)

args <- commandArgs(TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript script.R <input> <title> <output>\n")
}

input_file <- args[1]
plot_title <- args[2]
output_pdf <- args[3]

data <- read.delim(input_file, sep = "\t")

p <- ggplot(data, aes(x = Group, y = Diff*100, color = Group)) +
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.6, color = "grey") +
  geom_jitter(width = 0.2, size = 2) +
  scale_color_manual(values = c("HCC" = "#CD2626", "Control" = "grey50")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
	axis.ticks.x = element_line(color = "black", size = 0.2),
	axis.ticks.y = element_line(color = "black", size = 0.2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5)) +
  labs(x = "Group", y = plot_title, title = plot_title) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, comparisons = list(c("Control", "HCC")))

ggsave(output_pdf, plot = p, width = 3, height = 6)
