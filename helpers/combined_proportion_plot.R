library(cowplot)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

proportion_barplot <- function(seurat_obj, celltype_col, condition_col, condition_title = condition_col,
                               condition_cols = c("darkgreen", "darkred")) {
  plot_data <- seurat_obj@meta.data %>%
    group_by(sample, .data[[condition_col]], .data[[celltype_col]]) %>%
    summarize(cell_count = n(), .groups = 'drop') %>%
    group_by(sample, .data[[condition_col]]) %>%
    mutate(fraction = cell_count / sum(cell_count) * 100) %>%
    ungroup()
  
  custom_color <- condition_cols
  names(custom_color) <- levels(seurat_obj@meta.data[[condition_col]])
  custom_fill <- scales::alpha(condition_cols, 0.5)
  names(custom_fill) <- levels(seurat_obj@meta.data[[condition_col]])
  
  p <- ggplot(plot_data, aes_string(x = celltype_col, y = "fraction", fill = condition_col)) +
    scale_fill_manual(values = custom_fill, guide = 'none') +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.9)) +
    geom_point(aes_string(color = condition_col), position = position_dodge(width = 0.9), size = 1.6) +
    scale_color_manual(values = custom_color) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.background = element_blank(),
      strip.placement = "outside",
      panel.spacing.x = unit(1, "lines"),
      panel.spacing.y = unit(0.5, "lines"),
      panel.background = element_rect(fill = "white"),
      strip.text = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_text(angle = 0)
    ) +
    labs(
      x = NULL, y = "Fraction of cells (%)", fill = condition_title, color = condition_title
    )
  
  return(p)
}

grouped_proportion_barplot <- function(seurat_obj, celltype_col, cell_group_col, condition_col, condition_title = condition_col, condition_cols = c("darkgreen", "darkred"), panel_bg_cols = NULL) {
  plot_data <- seurat_obj@meta.data %>%
    group_by(sample, .data[[condition_col]], .data[[celltype_col]], .data[[cell_group_col]]) %>%
    summarize(cell_count = n(), .groups = 'drop') %>%
    group_by(sample, .data[[condition_col]]) %>%
    mutate(fraction = cell_count / sum(cell_count) * 100) %>%
    ungroup()
  
  bg_data <- unique(plot_data[, c(cell_group_col, celltype_col)])
  bg_data[[condition_col]] <- "dummy"
  bg_data$fraction <- 1
  
  custom_color <- condition_cols
  names(custom_color) <- levels(seurat_obj@meta.data[[condition_col]])
  custom_fill <- scales::alpha(condition_cols, 0.5)
  names(custom_fill) <- levels(seurat_obj@meta.data[[condition_col]])
  
  if (is.null(panel_bg_cols)) {
    unique_groups <- length(unique(seurat_obj@meta.data[[cell_group_col]]))
    panel_bg_cols <- brewer.pal(min(unique_groups, 9), "Pastel1")
    names(panel_bg_cols) <- unique(seurat_obj@meta.data[[cell_group_col]])
  } else {
    panel_bg_cols <- scales::alpha(panel_bg_cols, 0.2)
  }
  
  p <- ggplot(plot_data, aes_string(x = celltype_col, y = "fraction", fill = condition_col)) +
    geom_rect(data = bg_data, aes_string(fill = cell_group_col), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
    scale_fill_manual(values = c(custom_fill, panel_bg_cols), guide = 'none') +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.9)) +
    geom_point(aes_string(color = condition_col), position = position_dodge(width = 0.9), size = 1.6) +
    scale_color_manual(values = custom_color) +
    facet_grid(as.formula(paste("~", cell_group_col)), scales = "free_x", space = "free") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.background = element_blank(),
      strip.placement = "outside",
      panel.spacing.x = unit(1, "lines"),
      panel.spacing.y = unit(0.5, "lines"),
      panel.background = element_rect(fill = "white"),
      strip.text = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_text(angle = 0)
    ) +
    labs(
      x = NULL, y = "Fraction of cells (%)", fill = condition_title, color = condition_title
    )
  
  return(p)
}

permutation_test_plot <- function(sc_utils_obj, celltype_col, condition_col,
                                  condition_title = condition_col, condition_cols = c("darkgreen", "darkred"),
                                  FDR_threshold = 0.05, log2FD_threshold = log2(1.5)) {
  meta_data <- as.data.frame(sc_utils_obj@meta_data)
  permutation_data <- as.data.frame(sc_utils_obj@results$permutation)
  
  plot_data <- permutation_data %>%
    left_join(meta_data, by = c("clusters" = celltype_col)) %>%
    mutate(significance = ifelse(FDR < FDR_threshold & abs(obs_log2FD) > log2FD_threshold, 
                                 paste("FDR <", FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)), 
                                 "n.s."),
           significance = factor(significance, levels = c(paste("FDR <", FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)), "n.s.")))
  
  plot_data$clusters <- factor(plot_data$clusters, levels = levels(meta_data[[celltype_col]]))
  
  p <- ggplot(plot_data, aes_string(x = "clusters", y = "obs_log2FD")) + 
    geom_pointrange(aes_string(ymin = "boot_CI_2.5", ymax = "boot_CI_97.5", color = "significance")) + 
    theme_bw() + 
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.background = element_blank(),
      strip.placement = "outside",
      panel.spacing.x = unit(1, "lines"),
      panel.spacing.y = unit(0.5, "lines"),
      panel.background = element_rect(fill = "white"),
      strip.text = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_text(angle = 0)
    ) +
    geom_hline(yintercept = log2FD_threshold, lty = 2) + 
    geom_hline(yintercept = -log2FD_threshold, lty = 2) + 
    geom_hline(yintercept = 0) + 
    scale_color_manual(values = c("salmon", "grey")) +
    labs(y = "Observed log2 Fold Difference", x = "", color = "Significance")
  
  return(p)
}

grouped_permutation_test_plot <- function(sc_utils_obj, celltype_col, cell_group_col, condition_col, condition_title = condition_col, condition_cols = c("darkgreen", "darkred"), panel_bg_cols = NULL, FDR_threshold = 0.05, log2FD_threshold = log2(1.5)) {
  meta_data <- as.data.frame(sc_utils_obj@meta_data)
  permutation_data <- as.data.frame(sc_utils_obj@results$permutation)
  
  plot_data <- permutation_data %>%
    left_join(meta_data, by = c("clusters" = celltype_col)) %>%
    mutate(significance = ifelse(FDR < FDR_threshold & abs(obs_log2FD) > log2FD_threshold, 
                                 paste("FDR <", FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)), 
                                 "n.s."),
           significance = factor(significance, levels = c(paste("FDR <", FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)), "n.s.")))
  
  plot_data$clusters <- factor(plot_data$clusters, levels = levels(meta_data[[celltype_col]]))
  
  bg_data <- unique(plot_data[, c(cell_group_col, "clusters")])
  bg_data$ymin <- -Inf
  bg_data$ymax <- Inf
  bg_data$obs_log2FD <- 0
  
  if (is.null(panel_bg_cols)) {
    unique_groups <- length(unique(plot_data[[cell_group_col]]))
    panel_bg_cols <- brewer.pal(min(unique_groups, 9), "Pastel1")
    names(panel_bg_cols) <- unique(plot_data[[cell_group_col]])
  } else {
    panel_bg_cols <- scales::alpha(panel_bg_cols, 0.2)
  }
  
  p <- ggplot(plot_data, aes_string(x = "clusters", y = "obs_log2FD")) + 
    geom_rect(data = bg_data, aes_string(fill = cell_group_col), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
    scale_fill_manual(values = panel_bg_cols, guide = 'none') +
    geom_pointrange(aes_string(ymin = "boot_CI_2.5", ymax = "boot_CI_97.5", color = "significance")) + 
    theme_bw() + 
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.background = element_blank(),
      strip.placement = "outside",
      panel.spacing.x = unit(1, "lines"),
      panel.spacing.y = unit(0.5, "lines"),
      panel.background = element_rect(fill = "white"),
      strip.text = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_text(angle = 0)
    ) +
    geom_hline(yintercept = log2FD_threshold, lty = 2) + 
    geom_hline(yintercept = -log2FD_threshold, lty = 2) + 
    geom_hline(yintercept = 0) + 
    scale_color_manual(values = c("salmon", "grey")) +
    labs(y = "Observed log2 Fold Difference", x = "", color = "Significance") +
    facet_grid(as.formula(paste("~", cell_group_col)), scales = "free_x", space = "free")
  
  return(p)
}

combined_permutation_test_plot <- function(seurat_obj, prop_test,
                                           celltype_col, condition_col,
                                           condition_title = condition_col,
                                           condition_cols = c("darkgreen", "darkred"),
                                           FDR_threshold = 0.05,
                                           log2FD_threshold = log2(1.5)) {
  comparative_plot <- proportion_barplot(seurat_obj, celltype_col, condition_col, condition_title, condition_cols)
  permutation_plot <- permutation_test_plot(prop_test, celltype_col, condition_col, condition_title, condition_cols, FDR_threshold, log2FD_threshold)
  combined_plot <- plot_grid(comparative_plot, permutation_plot, ncol = 1, align = "v", axis = "lr")
  
  return(combined_plot)
}

grouped_combined_permutation_test_plot <- function(seurat_obj, prop_test, celltype_col, cell_group_col, condition_col, condition_title = condition_col, condition_cols = c("darkgreen", "darkred"), panel_bg_cols = NULL, FDR_threshold = 0.05, log2FD_threshold = log2(1.5)) {
  comparative_plot <- grouped_proportion_barplot(seurat_obj, celltype_col, cell_group_col, condition_col, condition_title, condition_cols, panel_bg_cols)
  permutation_plot <- grouped_permutation_test_plot(prop_test, celltype_col, cell_group_col, condition_col, condition_title, condition_cols, panel_bg_cols, FDR_threshold, log2FD_threshold)
  combined_plot <- plot_grid(comparative_plot, permutation_plot, ncol = 1, align = "v", axis = "lr")
  
  return(combined_plot)
}
