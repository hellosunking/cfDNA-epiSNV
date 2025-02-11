library(ggplot2)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript script.R <input> <title> <output_prefix>\n")
}

input_file <- args[1]
title_prefix <- args[2]
output <- args[3]

data <- read.table(input_file, header = TRUE, sep = "\t")

groups <- unique(data$Group)

output_pdf <- paste0(output, ".methylation.pdf")

pdf(output_pdf, width = 4, height = 6)

for (group in groups) {
  group_data <- data[data$Group == group, ]
  
  title_final <- paste0(group, ".", title_prefix)  

  group_data <- group_data[, c("Sample", "Group", "Wt_DNA", "Mut_DNA")]

  group_data_long <- melt(group_data, id.vars = c("Sample", "Group"), variable.name = "Type", value.name = "Value")
  
  group_data_long$Position <- ifelse(group_data_long$Type == "Wt_DNA", 1, 2)
  
  ref_values <- group_data_long$Value[group_data_long$Type == "Wt_DNA"]
  alt_values <- group_data_long$Value[group_data_long$Type == "Mut_DNA"]
  t_test_result <- t.test(ref_values, alt_values, paired = TRUE)
  
  p <- ggplot(group_data_long, aes(x = Position, y = Value * 100, group = Sample)) +
    geom_line(aes(group = Sample), colour = "grey70", size = 0.3) + 
    geom_point(aes(color = Type), size = 2, shape = 1, stroke = 1) + 
    scale_color_manual(values = c("Wt_DNA" = "grey50", "Mut_DNA" = "red")) +
    theme_minimal() +
    labs(x = "", y = "Percent(%)", title = title_final) +
    scale_x_continuous(breaks = c(1, 2), labels = c("Wt_DNA", "Mut_DNA"), limits = c(0.8, 2.2)) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 1),
      axis.ticks.x = element_line(color = "black", size = 0.2),
      axis.ticks.y = element_line(color = "black", size = 0.2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      plot.margin = margin(1, 1, 1, 1, "cm")
    ) +
    annotate("text", x = 1.5, y = max(group_data_long$Value * 100), 
             label = paste("p-value:", format(t_test_result$p.value, digits = 4)), 
             size = 4, hjust = 0.5)
  
  print(p)
}

dev.off()
