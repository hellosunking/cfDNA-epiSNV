library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript script.R <input> <title> <output>\n")
}

input_file <- args[1]
plot_title <- args[2]
output_pdf <- args[3]

data <- read.delim(input_file, sep = "\t")

data$Stage <- factor(data$Stage, levels = c("Control", "StageI", "StageII", "StageIII", "StageIV", "StageNA"))

colors <- c("Control" = "grey50", "StageI" = "deepskyblue", "StageII" = "darkorange4", "StageIII" = "darkorange", "StageIV" = "#CD2626", "StageNA" = "turquoise4")

p <- ggplot(data, aes(x = Stage, y = pred, color = Stage)) +
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.6, color = "grey") +
  geom_jitter(width = 0.2, size = 2) +
  scale_color_manual(values = colors) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
	axis.ticks.x = element_line(color = "black", size = 0.2),
	axis.ticks.y = element_line(color = "black", size = 0.2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "black"),
        plot.title = element_text(size = 14, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 90, hjust = 1)) +
  labs(x = "Stage", y = "Pred", title = plot_title) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, comparisons = list(c("Control", "StageI"), c("Control", "StageII"), c("Control", "StageIII"), c("Control", "StageIV"), c("Control", "StageNA")))

ggsave(output_pdf, plot = p, width = 12, height = 6)
