library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
plot_title <- args[2]
output_pdf <- args[3]

data <- read.delim(input_file, sep = "\t")

data$Group <- factor(data$Group, levels = c("Ctr", "HCC"))

p <- ggplot(data, aes(x = Group, y = sbs_contribution * 100, color = Group)) +
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.6, color = "grey") +
  geom_jitter(width = 0.2, size = 2) +
  #geom_quasirandom(groupOnX = TRUE, size = 1.5, dodge.width = 0.4, bandwidth = 0.1) +
  scale_color_manual(values = c("Ctr" = "grey50", "HCC_A" = "blue", "HCC_B" = "darkorange4", "HCC" = "#CD2626")) +
  theme_minimal(base_size = 14) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        plot.title = element_text(size = 14, hjust = 0.5)) +
  labs(x = "Group", y = "SBS_contribution(%)", title = plot_title) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, comparisons = list(c("Ctr", "HCC")))

ggsave(output_pdf, plot = p, width = 3, height = 6)
