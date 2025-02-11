library(ggplot2)
library(reshape2)

args <- commandArgs(TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript script.R <input> <title> <group> <output>\n")
}

input_file <- args[1]  
plot_title <- args[2]
group_name <- args[3]  
output_file <- args[4]  

data <- read.table(input_file, header = TRUE, sep = "\t")

# Filter data based on the group name
data <- subset(data, Group == group_name)

# Keep only columns 1, 3, and 4
data <- data[, c(1, 3, 4)]

# Melt the data to long format
data_long <- melt(data, id.vars = "Sample", variable.name = "Type", value.name = "Value")

# Add Position column based on Type
data_long$Position <- ifelse(data_long$Type == "Wt_DNA", 1, 2)

# Perform paired t-test on Wt_DNA and Mut_DNA values
ref_values <- data_long$Value[data_long$Type == "Wt_DNA"]
alt_values <- data_long$Value[data_long$Type == "Mut_DNA"]
t_test_result <- t.test(ref_values, alt_values, paired = TRUE)

# Create the plot
ggplot(data_long, aes(x = Position, y = Value * 100, group = Sample)) +
  geom_line(aes(group = Sample), colour = "grey70", size = 0.3) +  # Change linewidth to size
  geom_point(aes(color = Type), size = 2, shape = 1, stroke = 1) + 
  scale_color_manual(values = c("Wt_DNA" = "grey50", "Mut_DNA" = "red")) +
  theme_minimal() +
  labs(x = "", y = "Percent(%)", title = plot_title) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Wt_DNA", "Mut_DNA"), limits = c(0.8, 2.2)) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.ticks.x = element_line(color = "black", size = 0.2),
    axis.ticks.y = element_line(color = "black", size = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Change linewidth to size
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  annotate("text", x = 1.5, y = max(data_long$Value * 100), 
           label = paste("p-value:", format(t_test_result$p.value, digits = 4)), 
           size = 4, hjust = 0.5)

# Save the plot
ggsave(output_file, width = 4, height = 6)
