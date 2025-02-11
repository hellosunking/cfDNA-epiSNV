library(ggplot2)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript script.R <input_file> <plot_title> <output_file>\n")
}

input_file <- args[1]
plot_title <- args[2]
output_file <- args[3]

custom_title <- "Single Dataset" 

data <- read.table(input_file, header = TRUE, sep = "\t")

group_order <- colnames(data)[-1]

data_long <- pivot_longer(data, cols = -ID, names_to = "Group", values_to = "Value")

data_long$Group <- factor(data_long$Group, levels = group_order)

data_long$File <- custom_title

p <- ggplot(data_long, aes(x = Group, y = Value * 100)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, fill = "gray") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black")) +
  labs(title = plot_title, x = "SBS", y = "relative_contribution(%)") +
  ylim(0, 1.8) +
  facet_grid(File ~ ., scales = "free_y") + 
  theme(strip.placement = "outside",
        strip.text.y.right = element_text(angle = 90),
        strip.background = element_blank())

ggsave(output_file, plot = p, device = "pdf", width = 30, height = 15, units = "cm")
