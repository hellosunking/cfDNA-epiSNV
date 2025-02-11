library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

plot_contribution <- function(contribution,
                              signatures = NA,
                              index = NA,
                              coord_flip = FALSE,
                              mode = c("relative", "absolute"),
                              palette = NA) {
  # Match argument
  mode <- match.arg(mode)
  
  contribution <- contribution[, index, drop = FALSE]
  
  # Initialize variables to avoid R CMD check complaints
  Sample <- Contribution <- Signature <- Group <- NULL
  
  # Calculate signature contribution in absolute number of signatures if required
  if (mode == "absolute" && !is.na(signatures) && !is.null(signatures)) {
    total_signatures <- colSums(signatures)
    contribution <- sweep(contribution, 2, total_signatures, FUN = "*")
  }
  
  # Prepare the data frame
  tb <- contribution %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Signature") %>%
    tidyr::pivot_longer(-Signature, names_to = "Sample", values_to = "Contribution") %>%
    dplyr::mutate(
      Sample = factor(Sample, levels = unique(Sample)),
      Signature = factor(Signature, levels = unique(Signature)),
      Color = case_when(
        str_detect(Signature, "^A") ~ "Control",
        str_detect(Signature, "^SBS") ~ "a", 
        TRUE ~ "Other"                          
      )
    )
  
  # Setup the plot geometry and labels for the main plot
  bar_geom <- if (mode == "absolute") {
    geom_bar(stat = "identity", colour = "black")
  } else {
    geom_bar(position = "fill", stat = "identity", colour = "black")
  }
  
  y_lab <- if (mode == "absolute") "Absolute contribution \n (no. mutations)" else "Relative contribution"
  
  # Main plot
  plot_main <- ggplot(tb, aes(x = Sample, y = Contribution, fill = Color)) +
    bar_geom +
    labs(x = "", y = y_lab) +
    theme_bw() +
    theme(
      text = element_text(size = 12), # Increase text size
      axis.text.x = element_blank(),  # Hide x axis labels
      axis.text.y = element_text(size = 14),
      axis.ticks.x = element_line(color = "black", size = 0.2),
      axis.ticks.y = element_line(color = "black", size = 0.2),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.title = element_blank(),
      panel.border = element_blank(), # Remove plot border
      axis.ticks = element_blank(),   # Remove ticks
      axis.line = element_line(color = "black", linewidth = 0.5) # Draw x and y axis lines
    ) +
    scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, by = 0.2)) # y-axis breaks every 0.2
  
  # If custom palette is provided, apply it
  if (!is.na(palette) && !is.null(palette)) {
    plot_main <- plot_main + scale_fill_manual(values = palette)
  } else {
    # Default color mapping
    plot_main <- plot_main + scale_fill_manual(values = c("Control" = "#dae3f3", "a" = "#fbe5d6", "Other" = "#d3d3d3"))
  }
  
  # Return the two separate plots
  list(main_plot = plot_main)
}

args <- commandArgs(TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript script.R <input_file> <output_file>\n")
}

input.file <- args[1]
output.file <- args[2]

a <- read.table(input.file, head=T, row.names=1, sep="\t")

plots <- plot_contribution(a, coord_flip = FALSE, mode = "relative", index = c(1:72))

pdf(file= output.file )

print(plots$main_plot)

dev.off()
