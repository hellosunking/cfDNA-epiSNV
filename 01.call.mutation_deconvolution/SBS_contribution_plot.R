library(ggplot2)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
output_file <- args[2]
plot_title <- args[3]

data <- read.table(input_file, header = TRUE, sep = "\t")

group_order <- colnames(data)[-1]

data_long <- pivot_longer(data, cols = -ID, names_to = "Group", values_to = "Value")

data_long$Group <- factor(data_long$Group, levels = group_order)

p <- ggplot(data_long, aes(x = Group, y = Value*100)) +
  geom_boxplot(outlier.shape = NA, color = "gray") +
  geom_jitter(width = 0.2, color = "grey50", size = 0.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black")) +
  labs(title = plot_title, x = "SBS", y = "relative_contribution(%)")

ggsave(output_file, plot = p, device = "pdf", width = 30, height = 10, units = "cm")
