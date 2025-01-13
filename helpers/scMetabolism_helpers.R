library(scMetabolism)
library(wesanderson)
library(ggplot2)
library(patchwork)
library(grid)
library(cowplot)

# this modifies the original function to use umapharm and allow split.by

scMetabolismDimplot <- function(obj, pathway, split.by = NULL, ncol = NULL, combine = TRUE) {
  size <- 0.2
  #cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")

  # UMAP coordinates
  umap.loc <- obj@reductions$umapharm@cell.embeddings
  row.names(umap.loc) <- colnames(obj)

  # Metabolism signature expression
  signature_exp <- obj@assays$METABOLISM$score
  input.pathway <- pathway
  signature_ggplot <- data.frame(umap.loc, t(signature_exp[input.pathway, ]))

  # Define axis limits for consistency
  xlim <- range(signature_ggplot$umapharm_1)
  ylim <- range(signature_ggplot$umapharm_2)
  color_lim <- range(signature_ggplot[, 3])

  # Add splitting by group
  if (!is.null(split.by)) {
    split <- FetchData(object = obj, vars = split.by, clean = TRUE)[split.by]
    signature_ggplot <- cbind(signature_ggplot, split)
  }

  # Plotting
  pal <- wes_palette("Zissou1", 100, type = "continuous")

  # Create individual plots for each split
  plots <- lapply(unique(signature_ggplot[[split.by]]), function(group) {
    data_subset <- signature_ggplot[signature_ggplot[[split.by]] == group, ]
    ggplot(data = data_subset, aes(x = umapharm_1, y = umapharm_2, color = data_subset[, 3])) + 
      geom_point(size = size) + 
      scale_fill_gradientn(colours = pal, limits = color_lim) + 
      scale_color_gradientn(colours = pal, limits = color_lim) + 
      labs(color = input.pathway, title = group) + 
      xlab("umapharm_1") + 
      ylab("umapharm_2") + 
      xlim(xlim) + 
      ylim(ylim) + 
      theme(aspect.ratio = 1) + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            legend.position = "none")
  })

  # Generate the combined plot without legend
  combined_plot <- wrap_plots(plots, ncol = ncol %||% length(unique(signature_ggplot[[split.by]])))

  # Create the legend plot from the first plot
  legend_plot <- ggplot(signature_ggplot, aes(x = umapharm_1, y = umapharm_2, color = signature_ggplot[, 3])) +
    geom_point(size = size) +
    scale_fill_gradientn(colours = pal, limits = color_lim) + 
    scale_color_gradientn(colours = pal, limits = color_lim) + 
    labs(color = input.pathway) +
    theme(legend.position = "right")

  # Extract the legend
  tmp <- ggplot_gtable(ggplot_build(legend_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]

  # Combine the plot with the legend
  combined_plot_with_legend <- plot_grid(combined_plot, NULL, legend, ncol = 3, rel_widths = c(0.55, 0.02, 0.2))

  # Return the final plot
  return(combined_plot_with_legend)
}

