library(ggplot2)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)  
input_file <- args[1]  
plot_title <- args[2]  
output_file <- args[3]  

data <- read.table(input_file, header = TRUE, sep = "\t")

data_long <- melt(data, id.vars = "Sample", variable.name = "Type", value.name = "Value")

data_long$Position <- ifelse(data_long$Type == "Wt", 1, 2)

ref_values <- data_long$Value[data_long$Type == "Wt"]
alt_values <- data_long$Value[data_long$Type == "Mut"]
t_test_result <- t.test(ref_values, alt_values, paired = TRUE)

ggplot(data_long, aes(x = Position, y = Value*100, group = Sample)) +
  geom_line(aes(group = Sample), colour = "grey70", linewidth = 0.3) + 
  geom_point(aes(color = Type), size = 2, shape = 1, stroke = 1) + 
  scale_color_manual(values = c("Wt" = "grey50", "Mut" = "red")) +
  theme_minimal() +
  labs(x = "", y = "Abs 50 bp rate (%)", title = plot_title) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Wt", "Mut"), limits = c(0.8, 2.2)) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.ticks.x = element_line(color = "black", size = 0.2),
    axis.ticks.y = element_line(color = "black", size = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  annotate("text", x = 1.5, y = max(data_long$Value*100), 
           label = paste("p-value:", format(t_test_result$p.value, digits = 4)), 
           size = 4, hjust = 0.5)

ggsave(output_file, width = 4, height = 6, dpi = 300)
