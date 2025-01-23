library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggalluvial)
library(colorspace)
library(cowplot)
library(DESeq2)
library(networkD3)
library(webshot2)
library(htmlwidgets)
library(png)
library(grid)
library(Seurat)
library(Platypus)
library(gridExtra)
library(forcats)
library(broom)
library(RColorBrewer)

cells_clonotype_size <- function(vgm) {
  vgm <- transfer_vdj_to_gex(vgm, "clonotype_id")
  clonotypes_table <- table(vgm[["GEX"]]$clonotype_id)
  vgm[["GEX"]]$clonotype_size <-
    clonotypes_table[as.character(vgm[["GEX"]]$clonotype_id)]
  vgm[["GEX"]]$clonotype_size[is.na(vgm[["GEX"]]$clonotype_size)] <- 1
  vgm <- transfer_gex_to_vdj(vgm, "clonotype_size")
  return(vgm)
}

clonotype_freq_plot <- function(vdj, relative=T, by = "sample",
                                bins = c(1, 1, 5, 20, Inf),
                                labels = c("1", "2-5", "6-20", ">20"),
                                bin_cols = c("#d73027", "#fc8d59",
                                             "#fee090", "#313695")) {
  clonotype_freq <- vdj %>% group_by(.data[[by]], clonotype_id)

  clonotype_freq <- clonotype_freq %>%
    summarise(count = n(), .groups = 'drop') %>%
    arrange(.data[[by]], desc(count))

  clonotype_freq <- clonotype_freq %>%
    group_by(.data[[by]])

  if (relative) {
    clonotype_freq <- clonotype_freq %>%
      mutate(cumulative_proportion = cumsum(count) / sum(count))
  }

  clonotype_freq <- clonotype_freq %>%
    mutate(bin = cut(row_number(), breaks = bins, labels = labels,
                     right = FALSE))

  plot_data <- clonotype_freq %>% group_by(.data[[by]], bin)

  plot_data <- plot_data %>%
    summarise(
      total_count = sum(count),
      proportion = sum(count) / sum(clonotype_freq$count), .groups = 'drop') %>%
    arrange(.data[[by]], bin)

  if (relative) {
    plt <- ggplot(plot_data, aes(x = .data[[by]],
                                 y = proportion, fill = bin)) +
      geom_bar(stat = "identity", position = "fill") +
      labs(y = "Clonal Proportion")
  } else {
    plt <- ggplot(plot_data, aes(x = .data[[by]],
                                 y = total_count, fill = bin)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(y = "Clonal Count")
  }

  plt + scale_fill_manual(values = bin_cols) +
    labs(x = "Samples", fill = "Top N clones") +
    theme_minimal() + theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

cmp_clonotype_size_plot <- function(vgm, relative, ident, group, bin_cols,
                                    sel_celltype = NA,
                                    celltype_col = "celltype.abbrev",
                                    compute_pval = F,
                                    hide_legend = F) {

  if (is.na(sel_celltype)) {
    selected_data <- vgm$GEX@meta.data
  } else {
    selected_data <-
      vgm$GEX@meta.data[vgm$GEX@meta.data[[celltype_col]] == sel_celltype, ]
  }

  plot_data <- selected_data %>%
    group_by(.data[[ident]], .data[[group]], clonotype_size_bin) %>%
    summarise(count = n(), .groups = 'drop') %>%
    arrange(.data[[ident]], .data[[group]], clonotype_size_bin)

  if (relative) {
    plot_data <- plot_data %>%
      group_by(.data[[ident]], .data[[group]]) %>%
      mutate(proportion = count / sum(count))
    plt <- ggplot(plot_data, aes(x = .data[[group]], y = proportion,
                                 fill = clonotype_size_bin)) +
      geom_bar(stat = "identity", position = "fill") +
      labs(y = "Clonotype Size (relative)") +
      scale_y_continuous(expand = c(0, 0.02), breaks = seq(0, 1.45, 0.25))

  } else {
    plt <- ggplot(plot_data, aes(x = .data[[group]], y = count,
                                 fill = clonotype_size_bin)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(y = "Clonotype Size (count)")
  }

  if (compute_pval) {
    test_results <- plot_data %>%
      group_by(.data[[ident]]) %>%
      summarise(test_result =
                list(broom::tidy(chisq.test(matrix(count, ncol = 2))))) %>%
      unnest(test_result) %>%
      select(.data[[ident]], statistic, p.value) %>%
      mutate(adj_p_value = p.adjust(p.value, method = "BH"))
  }

  plt = plt +
    scale_fill_manual(values = bin_cols) +
    labs(fill = "Clonotype Size") +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 25)) +
    theme(axis.text.y = element_text(size = 16)) +
    theme(plot.title = element_text(size=30, hjust=0.5)) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())

  if (is.na(sel_celltype)) {
    plt = plt + facet_wrap(~.data[[ident]], ncol = 2)
  } else {
    plt = plt + ggtitle(sel_celltype)
  }

  if (compute_pval) {
    if (relative) {
      plt = plt +
        geom_text(data = test_results,
            aes(x = 1.5, y = 1.07,
                label = ifelse(adj_p_value > 0.05,
                           paste0("n.s. (", sprintf("%.3f", adj_p_value), ")"),
                           ifelse(adj_p_value > 0.001,
                                  sprintf("%.3f", adj_p_value),
                                  ifelse(adj_p_value < 1e-323,
                                         "< 1.000 e-323",
                                         sprintf("%.3e", adj_p_value))))),
            inherit.aes = FALSE, size = 8) +
        geom_segment(data = test_results,
                     aes(x = 1, xend = 2,
                         y = 1.03, yend = 1.03),
                     inherit.aes = FALSE, size = 2, color="black")

    } else {
      plt = plt +
        geom_text(data = test_results,
            aes(x = 1.5, y = max(plot_data$count)*1.1,
                label = ifelse(adj_p_value > 0.05,
                           paste0("n.s. (", sprintf("%.3f", adj_p_value), ")"),
                           ifelse(adj_p_value > 0.001,
                                  sprintf("%.3f", adj_p_value),
                                  ifelse(adj_p_value < 1e-323,
                                         "< 1.000 e-323",
                                         sprintf("%.3e", adj_p_value))))),
            inherit.aes = FALSE, size = 8) +
        geom_segment(data = test_results,
                     aes(x = 1, xend = 2,
                         y = max(plot_data$count)*1.05,
                         yend = max(plot_data$count)*1.05),
                     inherit.aes = FALSE, size = 2, color="black")
    }
  }

  if (hide_legend) {
    plt = plt + theme(legend.position = "none")
  }

  plt
}

celltype_cmp_clonotype_size_plot <- function(vgm, sel_celltype,
                                             ident_group_table,
                                             celltype_col = "celltype",
                                             ident_col,
                                             group_col, group1, group2,
                                             bin_cols) {
  selected_data <-
    vgm$GEX@meta.data[vgm$GEX@meta.data[[celltype_col]] == sel_celltype, ]

  plot_data <- selected_data %>%
    group_by(.data[[ident_col]], clonotype_id,
             .data[[celltype_col]], .data[[group_col]]) %>%
    summarise(count = n(), .groups = 'drop') %>%
    arrange(.data[[ident_col]], desc(count))

  plot_data <- plot_data %>%
    mutate(bin =
           factor(vgm$GEX$clonotype_size[match(clonotype_id,
                                                   vgm$GEX$clonotype_id)])) %>%
   filter(!is.na(bin))

  plot_data <- plot_data %>%
    group_by(.data[[ident_col]], .data[[group_col]], bin) %>%
    summarise(total_count = sum(count), .groups = 'drop') %>%
    mutate(proportion = total_count / sum(total_count)) %>%
    arrange(.data[[ident_col]], .data[[group_col]], bin)

  # Rename the clinical_condition column in sample_cond_table to avoid conflict
  ident_group_table$group = ident_group_table[[group_col]]
  ident_group_table[[group_col]] = NULL

  # Join with sample condition table to ensure correct sorting
  plot_data <- plot_data %>%
    left_join(ident_group_table, by = ident_col) %>%
    arrange(factor(group, levels = c(group1, group2)), .data[[ident_col]])

  # Add a blank bar for separation
  blank_row <- tibble(
    sample = " ",
    group = NA,
    bin = factor("1", levels = levels(plot_data$bin)),
    total_count = 0,
    proportion = 0,
    sample_condition = NA
  )

  plot_data <- bind_rows(
    plot_data %>% filter(group == group1),
    blank_row,
    plot_data %>% filter(group == group2)
  )

  plot_data <- plot_data %>% filter(!is.na(.data[[ident_col]]))

  n_samples <- n_distinct(plot_data[[ident_col]])
  n_group1_samples <- plot_data %>%
    filter(.data[[group_col]] == group1) %>%
    summarise(count = n_distinct(.data[[ident_col]])) %>%
    pull(count)
  n_group2_samples <- plot_data %>%
    filter(.data[[group_col]] == group2) %>%
    summarise(count = n_distinct(.data[[ident_col]])) %>%
    pull(count)

  ggplot(plot_data, aes(x = factor(.data[[ident_col]],
                        levels = unique(.data[[ident_col]])),
                        y = proportion, fill = bin)) +
    geom_bar(stat = "identity", position = "fill") +
    ggtitle(sel_celltype) +
    scale_fill_manual(values = bin_cols) +
    labs(y = "Clonal Proportion", fill = "Clonotype Size") +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_segment(aes(x = 1, xend = n_group1_samples,
                     y = -0.05, yend = -0.05),
                 size = 1.5, color = "black") +
    geom_segment(aes(x = n_group1_samples + 2, xend = n_samples,
                     y = -0.05, yend = -0.05),
                 size = 1.5, color = "black") +
    annotate("text", x = 0.5 + n_group1_samples/2,
             y = -0.1, label = group1, size = 5, hjust = 0.5) +
    annotate("text", x = 1.5 + n_group1_samples + n_group2_samples/2 ,
             y = -0.1, label = group2, size = 5, hjust = 0.5)
}

palette_DJC_locus =  c(brewer.pal(9, "Blues")[5], brewer.pal(9, "Reds")[4])
palette_J_gene = c(brewer.pal(9, "Set1"),brewer.pal(8, "Set2"),
                   brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"))

cols_vdj = c(brewer.pal(9, "Set1"),brewer.pal(8, "Set2"),
             brewer.pal(12, "Set3")[3:9])
cols_vj = c(brewer.pal(9, "Set1"),brewer.pal(8, "Set2"),
            brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"))

transfer_gex_to_vdj <- function(vgm, cols) {
  for (col in cols) {
    vgm[["VDJ"]][[col]] <- NULL
  }
  vgm[["VDJ"]] <- merge(vgm[["VDJ"]],
                        vgm[["GEX"]]@meta.data[, c("barcode", cols)],
                        all.x = T, by="barcode", sort=F, no.dups=F)
  vgm
}

transfer_vdj_to_gex <- function(vgm, cols) {
  for (col in cols) {
    vgm[["GEX"]]@meta.data[[col]] <- NULL
  }
  newmeta <- merge(vgm[["GEX"]]@meta.data,
                   vgm[["VDJ"]][, c("barcode", cols)],
                   all.x = T, by="barcode", sort=F, no.dups=F)
  rownames(newmeta) <- newmeta$barcode
  vgm[["GEX"]]@meta.data <- newmeta
  vgm
}

vgm_drop_levels <- function(vgm) {
  for (col in colnames(vgm[["GEX"]]@meta.data)) {
    if (is.factor(vgm[["GEX"]]@meta.data[[col]])) {
      vgm[["GEX"]]@meta.data[[col]] <- droplevels(vgm[["GEX"]]@meta.data[[col]])
    }
  }
  for (col in colnames(vgm[["VDJ"]])) {
    if (is.factor(vgm[["VDJ"]][[col]])) {
      vgm[["VDJ"]][[col]] <- droplevels(vgm[["VDJ"]][[col]])
    }
  }
  return(vgm)
}

vgm_filter_vj_too_long <- function(vgm, max_length) {
  vgm$VDJ <- vgm$VDJ[nchar(vgm$VDJ$VJ_cdr3s_aa) <= max_length, ]
  vgm[["GEX"]]$VDJ_available <-
    vgm[["GEX"]]$barcode %in% vgm[["VDJ"]]$barcode
  vgm <- vgm_drop_levels(vgm)
  return(vgm)
}

vgm_filter_unavailable <- function(vgm) {
  vgm[["VDJ"]] <- subset(vgm[["VDJ"]], GEX_available == T)
  vgm[["GEX"]] <- subset(vgm[["GEX"]], VDJ_available == T)
  vgm <- vgm_drop_levels(vgm)
  return(vgm)
}

require_at_least_one_H_and_one_L <- function(vgm) {
  vgm[["VDJ"]] <-
    subset(vgm[["VDJ"]], Nr_of_VDJ_chains >= 1 & Nr_of_VJ_chains >= 1)
  vgm[["GEX"]]$VDJ_available <-
    vgm[["GEX"]]$barcode %in% vgm[["VDJ"]]$barcode
  vgm[["GEX"]] <- subset(vgm[["GEX"]], VDJ_available == T)
  vgm <- vgm_drop_levels(vgm)
  return(vgm)
}

require_one_H_and_not_more_than_one_L <- function(vgm) {
  vgm[["VDJ"]] <-
    subset(vgm[["VDJ"]], Nr_of_VDJ_chains == 1 & Nr_of_VJ_chains <= 1)
  vgm[["GEX"]]$VDJ_available <-
    vgm[["GEX"]]$barcode %in% vgm[["VDJ"]]$barcode
  vgm[["GEX"]] <- subset(vgm[["GEX"]], VDJ_available == T)
  vgm <- vgm_drop_levels(vgm)
  return(vgm)
}

vgm_filter_by_celltype <- function(vgm, keep, column) {
  vgm[["VDJ"]] <- subset(vgm[["VDJ"]], eval(parse(text = column)) %in% keep)
  vgm[["GEX"]]$VDJ_available <-
    vgm[["GEX"]]$barcode %in% vgm[["VDJ"]]$barcode
  vgm[["GEX"]] <- subset(vgm[["GEX"]], VDJ_available == T)
  vgm <- vgm_drop_levels(vgm)
  return(vgm)
}

vgm_filter_unavailable_vj <- function(vgm) {
  vgm[["VDJ"]] <-
    subset(vgm[["VDJ"]], Nr_of_VJ_chains == 1)
  vgm[["GEX"]]$VDJ_available <-
    vgm[["GEX"]]$barcode %in% vgm[["VDJ"]]$barcode
  vgm[["GEX"]] <- subset(vgm[["GEX"]], VDJ_available == T)
  vgm <- vgm_drop_levels(vgm)
  return(vgm)
}

#
# This is a modified version of Platypus::GEX_phenotype.
# I modified the line before the return to:
#     Seurat::Idents(seurat.object) <- "previous.ident"
# since it gave an error, otherwise.
#
mod_GEX_phenotype <- function (seurat.object, cell.state.names,
                               cell.state.markers, default)
{
    if (missing(default))
        default <- TRUE
    Cap <- function(x) {
        temp <- c()
        for (i in 1:length(x)) {
            s <- strsplit(x, ";")[[i]]
            temp[i] <- paste(toupper(substring(s, 1, 1)), tolower(substring(s,
                2)), sep = "", collapse = ";")
        }
        return(temp)
    }
    is.hum <- any(useful::find.case(rownames(seurat.object),
        case = "upper"))
    if (missing(cell.state.markers) & default == T) {
        cell.state.markers <- c("CD4+;CD44-", "CD4+;IL7R+;CD44+",
            "CD4+;CD44+;IL7R-;IFNG+", "CD8A+;TCF7+;CD44-",
            "CD8A+;CX3CR1+;IL7R-", "CD8A+;IL7R+;CD44+", "PDCD1+;CD8A+",
            "CD19+;CD27-;CD38-", "FAS+;CD19+", "SDC1+", "CD38+;FAS-")
    }
    if (missing(cell.state.names) & default == T) {
        cell.state.names <- c("NaiveCd4", "MemoryCd4", "ActivatedCd4",
            "NaiveCd8", "EffectorCd8", "MemoryCd8", "ExhaustedCd8",
            "NaiveBcell", "GerminalcenterBcell", "Plasmacell",
            "MemoryBcell")
    }
    if (is.hum == F && default == T) {
        if (is.hum == F) {
            cell.state.markers <- Cap(cell.state.markers)
        }
        if (is.hum == T && default == T) {
            cell.state.markers <- toupper(cell.state.markers)
        }
    }
    cell.state.markers <- gsub(pattern = ";", replacement = "&",
        cell.state.markers)
    cell.state.markers <- gsub(pattern = "\\+", replacement = ">0",
        cell.state.markers)
    cell.state.markers <- gsub(pattern = "-", replacement = "==0",
        cell.state.markers)
    seurat.object[["previous.ident"]] <- Seurat::Idents(object = seurat.object)
    Seurat::Idents(seurat.object) <- "Unclassified"
    cmd <- c()
    for (i in 1:length(cell.state.names)) {
        cmd[i] <- paste0(cell.state.names[i],
                         "<-Seurat::WhichCells(seurat.object, ",
                         "slot = 'counts', expression =",
                         cell.state.markers[i], ")")
        is.exist <- tryCatch(expr = length(eval(parse(text = cmd[i]))),
            error = function(x) {
                x <- F
                return(x)
            })
        if (is.exist != F) {
            Seurat::Idents(object = seurat.object,
                           cells = eval(parse(text = cell.state.names[i]))) <-
                             cell.state.names[i]
        }
    }
    seurat.object[["cell.state"]] <- Seurat::Idents(object = seurat.object)
    Seurat::Idents(seurat.object) <- "previous.ident"
    return(seurat.object)
}

VDJ_overlap_sorted_heatmap <- function(VDJ, feature.columns, grouping.column,
                                       group.by, pvalues.label.size = 4,
                                       axis.label.size = 4)
{
    group <- NULL
    overlap <- NULL
    overlap_lab <- NULL
    platypus.version <- "v3"

    if (!"barcode" %in% names(VDJ)) {
        warning("'barcode' column must be present in input dataframe")
    }

    to_remove <- c()
    for (n in 1:nrow(VDJ)) {
        if ("" %in% VDJ[n, c(feature.columns)]) {
            to_remove <- c(to_remove, n)
        }
    }
    if (length(to_remove) > 0) {
        VDJ <- VDJ[-to_remove, ]
    }

    grouping <- data.frame(group = VDJ[, grouping.column])
    if (NA %in% grouping$group) 
        stop(paste0("NA values in grouping columns. ",
             "Please choose another column or replace NA values"))

    grouping$group <- factor(grouping$group, levels =
                             levels(VDJ[, grouping.column]))

    if (length(feature.columns) > 1) {
        grouping$pasted <-
          do.call(paste, c(VDJ[, c(feature.columns)], sep = "/"))
    } else {
        grouping$pasted <- VDJ[, c(feature.columns)]
    }

    # Order by group.by column, then by factor levels
    if (!missing(group.by)) {
        grouping$group.by <- VDJ[, group.by]
        grouping <- grouping[order(grouping$group.by, grouping$group), ]
    }

    sample.names <- unique(grouping[, 1])
    df.list <- list()
    for (i in 1:length(unique(grouping[, 1]))) {
        df.list[[i]] <-
          unique(subset(grouping, grouping[, 1] ==
                        unique(grouping[, 1])[i])[, 2])
    }
    names(df.list) <- sample.names
    print(sample.names)

    if (length(sample.names) > 2) {
        combs <-
          as.data.frame(t(utils::combn(as.character(sample.names),
                                       m = 2, simplify = TRUE)))
        combs[, 1] <- factor(combs[, 1], levels = sample.names)
        combs[, 2] <- factor(combs[, 2], levels = sample.names)
    } else {
        combs <- data.frame(sample.names[1], sample.names[2])
    }

    combs$overlap <- NA
    combs$items.overlapping <- NA
    ov_temp_list <- list()
    for (i in 1:nrow(combs)) {
        if (all(is.na(df.list[[which(names(df.list) == combs[i, 1])]])) == F &
            all(is.na(df.list[[which(names(df.list) == combs[i, 2])]])) == F) {
            combs$overlap[i] <-
              sum(df.list[[which(names(df.list) == combs[i, 1])]] %in%
                  df.list[[which(names(df.list) == combs[i, 2])]])
            ov_temp <-
              df.list[[which(names(df.list) == combs[i, 1])]][
                which(df.list[[which(names(df.list) == combs[i, 1])]] %in%
                      df.list[[which(names(df.list) == combs[i, 2])]])]
            combs$items.overlapping[i] <- paste0(ov_temp, collapse = ";")
        } else {
            combs$overlap[i] <- NA
            ov_temp <- ""
        }
    }

    combs$overlap_lab <- round(combs$overlap, 3)
    combs$overlap_lab <- as.character(combs$overlap_lab)
    combs$overlap_lab[is.na(combs$overlap)] <- "NA"

    plot_out <- ggplot2::ggplot(combs, ggplot2::aes(x = combs[, 1],
                                                    y = combs[, 2],
                                                    fill = overlap)) +
        ggplot2::geom_tile() +
        ggplot2::geom_text(ggplot2::aes(label = overlap_lab),
                           size = pvalues.label.size) +
        ggplot2::scale_fill_gradient2(low = "navy", mid = "white", high = "red",
                                      limits = range(combs$overlap)) +
        ggplot2::theme(
            panel.background = ggplot2::element_blank(),
            axis.text = ggplot2::element_text(size = 30),
            axis.line.x = ggplot2::element_blank(),
            axis.line.y = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            text = ggplot2::element_text(size = 30),
            legend.key = ggplot2::element_rect(colour = "white"),
            legend.position = "none",
            plot.title = ggplot2::element_text(hjust = 0.5, size = 25),
            plot.subtitle = ggplot2::element_text(size = 15),
            axis.text.x = ggplot2::element_text(angle = 60, vjust = 1,
                                           hjust = 1, size = axis.label.size),
            axis.text.y = ggplot2::element_text(size = axis.label.size)
        ) +
        ggplot2::labs(
            title = "",
            x = "",
            y = "",
            subtitle = paste0("Overlap features: ",
                              paste0(feature.columns, collapse = " ; ")),
            fill = ""
        ) +
        ggplot2::scale_y_discrete(limits = rev)

    return(list(plot_out, combs))
}

most_freq_cdr3s_aa_plot <- function(vdj, use.VJ = FALSE, group_col =
                                   "clinical_condition", color_by = "sample",
                                   cols = colorRampPalette(brewer.pal(9,
                                     "Set1"))(length(unique(vdj[[color_by]]))),
                                   color_lbl = color_by,
                                   title = "Most frequent CD3 sequences",
                                   sel_celltype = NA, celltype_col = "celltype",
                                   top_n = 20, short_title=F) {
  if (use.VJ) {
    aa = "VJ_cdr3s_aa"
  } else {
    aa = "VDJ_cdr3s_aa"
  }
  if (!short_title) {
    if (use.VJ) {
      title = paste0(title, " (VJ)")
    } else {
      title = paste0(title, " (VDJ)")
    }
  }

  if (!is.na(sel_celltype)) {
    vdj <- vdj %>%
      filter(.data[[celltype_col]] == sel_celltype)
    if (nrow(vdj) == 0) {
      stop(paste("No data found for cell type:", sel_celltype))
    }
    if (title == "") {
      title = sel_celltype
    } else if (title == " (VJ)" | title == " (VDJ)") {
        title = paste0(sel_celltype, title)
    } else {
      title = paste0(sel_celltype, ": ", title)
    }
  }

  vdj <- vdj %>% filter(.data[[aa]] != "")

  # Add zero-width space to sequences in the second condition
  second_condition <- unique(vdj[[group_col]])[2]
  vdj <- vdj %>%
    mutate(!!aa := ifelse(.data[[group_col]] == second_condition,
                          paste0(.data[[aa]], "\u200B"),
                          .data[[aa]]))

  total_sequences <- vdj %>%
    group_by(.data[[group_col]]) %>%
    summarise(total_count = n(), .groups = "drop")

  # Merge total count back into the original data
  cdr3_usage <- vdj %>%
    group_by(.data[[group_col]], .data[[aa]], .data[[color_by]]) %>%
    summarise(count = n(),
              celltypes = paste(unique(.data[[celltype_col]]), collapse = ";"),
              .groups = "drop") %>%
    left_join(total_sequences, by = group_col) %>%
    mutate(mean_usage = count / total_count)

  # Summarize total usage for each sequence across all groups
  total_usage <- cdr3_usage %>%
    group_by(.data[[group_col]], .data[[aa]]) %>%
    summarise(seq_count = sum(count), .groups = "drop") %>%
    left_join(total_sequences, by = group_col) %>%
    mutate(seq_mean_usage = seq_count / total_count) %>%
    arrange(.data[[group_col]], desc(seq_mean_usage))

  # Get the top n sequences by total mean usage for each group
  top_sequences <- total_usage %>%
    group_by(.data[[group_col]]) %>%
    slice_head(n = top_n) %>%
    ungroup()

  # Filter cdr3_usage to keep only top sequences for each group
  cdr3_usage_top <- cdr3_usage %>%
    semi_join(top_sequences, by = c(group_col, aa))

  cdr3_usage_top$total_count <- NULL

  cdr3_usage_top <- cdr3_usage_top %>%
    left_join(total_usage, by = c(group_col, aa)) %>%
    arrange(.data[[group_col]], desc(seq_mean_usage), .data[[aa]])

  cdr3_usage_top <- cdr3_usage_top %>%
    group_by(.data[[group_col]]) %>%
    mutate(aa_reordered = fct_reorder(.data[[aa]], -seq_mean_usage)) %>%
    ungroup()

  plt <- ggplot(cdr3_usage_top, aes(x = aa_reordered, y = count,
                                    fill = .data[[color_by]])) +
    geom_bar(stat = "identity") +
    facet_wrap(~ .data[[group_col]], scales = "free_x", ncol = 2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.margin = margin(0,0,-100,0),
          axis.title.x = element_blank()) +
    labs(x = "CDR3 Sequences", y = "Count",
         title = title, fill = color_lbl) +
    scale_fill_manual(values = cols)

  retval = list()
  retval[["cdr3_usage_top"]] = cdr3_usage_top
  retval[["plot"]] = plt
  return(retval)
}

gene_usage_alluvial <- function(vdj,
                                left, right,
                                left_lbl = left, right_lbl = right,
                                left_lbl_pfx = "", right_lbl_pfx = "",
                                grp_col = "group_id",
                                grp1 = unique(vdj[[grp_col]])[[1]],
                                grp2 = unique(vdj[[grp_col]])[[2]],
                                sel_celltype = NA, celltype_col = "celltype",
                                min_clonotype_frequency = NA,
                                clonotype_frequency_col = "clonotype_frequency",
                                signif_pairs = NULL,
                                cols = NULL,
                                stratum_label_size = 2.4,
                                stratum_label_size_prop_factor = 12.5,
                                stratum_label_threshold = 0.01,
                                group_title_size = 14,
                                main_title_size = 18,
                                y_axis_title_size = 10,
                                y_axis_text_size = 8,
                                x_axis_text_size = 12,
                                alluvium_width = 0.25,
                                stratum_width = 0.25,
                                not_signif_col = "#E0E0E0") {

  if (!is.na(sel_celltype)) {
    vdj <- vdj %>%
      filter(.data[[celltype_col]] == sel_celltype)
    if (nrow(vdj) == 0) {
      stop(paste("No data found for cell type:", sel_celltype))
    }
  }

  if (!is.na(min_clonotype_frequency)) {
    vdj <- vdj %>%
      filter(.data[[clonotype_frequency_col]] >= min_clonotype_frequency)
    if (nrow(vdj) == 0) {
      stop(paste("No data found with clonotype size >=",
                 min_clonotype_frequency))
    }
  }

  alluvial_data <- vdj %>%
    dplyr::select(.data[[left]], .data[[right]], .data[[grp_col]]) %>%
    group_by(.data[[left]], .data[[right]], .data[[grp_col]]) %>%
    summarise(freq = n()) %>%
    ungroup()

  if (is.null(cols)) {
    # palette with no grays
    cols <- c(
              brewer.pal(9, "Set1")[c(1:8)],
              brewer.pal(8, "Set2")[c(1:7)],
              brewer.pal(12, "Set3")[c(1:8)],
              brewer.pal(12, "Set3")[c(10:12)],
              brewer.pal(8, "Dark2")[1:7]
    )
  }
  if (length(cols) < length(unique(vdj[[left]]))) {
    cols = colorRampPalette(cols)(length(unique(vdj[[left]])))
  }

  color_mapping <- setNames(cols, unique(alluvial_data[[left]]))

  if (!is.null(signif_pairs)) {
    alluvial_data <- alluvial_data %>%
      mutate(pair = paste(.data[[left]], .data[[right]], sep = '+'),
             signif_status = ifelse(pair %in% signif_pairs,
                                    "Significant", "Not Significant"),
             fill_color = ifelse(signif_status == "Significant",
                                 color_mapping[.data[[left]]],
                                 not_signif_col))
  } else {
    alluvial_data <- alluvial_data %>%
      mutate(fill_color = color_mapping[.data[[left]]])
  }

  if (left_lbl_pfx != "") {
    alluvial_data$left_plot <-
      sub(paste0("^", left_lbl_pfx), "", alluvial_data[[left]])
  } else {
    alluvial_data$left_plot <- alluvial_data[[left]]
  }
  if (right_lbl_pfx != "") {
    alluvial_data$right_plot <-
      sub(paste0("^", right_lbl_pfx), "", alluvial_data[[right]])
  } else {
    alluvial_data$right_plot <- alluvial_data[[right]]
  }
  alluvial_data_grp1 <- alluvial_data %>% filter(.data[[grp_col]] == grp1)
  alluvial_data_grp2 <- alluvial_data %>% filter(.data[[grp_col]] == grp2)

  gene_usage_alluvial_subplot <- function(data, title, colors) {
    ggplot(data, aes(axis1 = left_plot, axis2 = right_plot, y = freq)) +
      geom_alluvium(aes(fill = fill_color), width = alluvium_width) +
      geom_stratum(width = stratum_width) +
      geom_text(stat = "stratum",
                aes(label = after_stat(stratum),
                    size = ifelse(after_stat(prop) < stratum_label_threshold, 0,
                                  stratum_label_size *
                                    (1 + prop *
                                     stratum_label_size_prop_factor)))) +
      scale_size_identity() +
      scale_x_discrete(limits = c(left_lbl, right_lbl), expand = c(0, 0)) +
      scale_fill_identity() +
      theme_minimal() +
      theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = y_axis_title_size),
        axis.text.y = element_text(size = y_axis_text_size),
        axis.text.x = element_text(size = x_axis_text_size),
        plot.title = element_text(hjust = 0.5, face = "bold",
                                  size = group_title_size),
      ) +
      labs(title = title, y = "Count") +
      theme(legend.position = "none")
  }

  plot_grp1 <- gene_usage_alluvial_subplot(alluvial_data_grp1, grp1)
  plot_grp2 <- gene_usage_alluvial_subplot(alluvial_data_grp2, grp2)
  plt <- plot_grid(plot_grp1, plot_grp2, ncol = 2)

  if (!is.na(sel_celltype)) {
    plt <- plot_grid(
      ggdraw() + draw_label(sel_celltype, fontface = 'bold',
                            size = main_title_size, hjust = 0.5),
      plt, ncol = 1, rel_heights = c(0.1, 1)
    )
  }

  return(plt)
}

gene_pairs_usage_deseq <- function(vdj, ident_grp_table, ident_col, locus1,
                                   locus2, grp_col, grp1, grp2) {
  vdj_frequencies <- vdj %>%
    group_by(.data[[locus1]], .data[[locus2]], .data[[ident_col]]) %>%
    summarise(freq = n()) %>%
    ungroup()
  count_matrix <- vdj_frequencies %>% spread(key = .data[[ident_col]],
                                             value = freq, fill = 0)
  count_matrix[is.na(count_matrix)] <- 0
  counts <- as.matrix(count_matrix %>%
                        select(-.data[[locus1]], -.data[[locus2]]))
  rownames(counts) <-
    paste(count_matrix[[locus1]], count_matrix[[locus2]], sep = "_")
  coldata <- as.data.frame(ident_grp_table)
  coldata <- coldata[match(colnames(counts), coldata[[ident_col]]), ]
  print(coldata)
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                                design = as.formula(paste("~", grp_col)))
  dds <- DESeq(dds)
  res <- results(dds)
  res[[locus1]] <- sapply(strsplit(rownames(res), "_"), `[`, 1)
  res[[locus2]] <- sapply(strsplit(rownames(res), "_"), `[`, 2)
  res$status <- ifelse(res$padj < 0.05 & res$log2FoldChange > 0,
                       paste0("Enriched in ", grp2),
                       ifelse(res$padj < 0.05 & res$log2FoldChange < 0,
                              paste0("Depleted in ", grp2), "Not Significant"))
  res
}

single_genes_usage_deseq <- function(vdj, features, ident_grp_table,
                                     ident_col = "sample_id",
                                     grp_col = "group_id",
                                     grp1 = unique(vdj[[grp_col]])[[1]],
                                     grp2 = unique(vdj[[grp_col]])[[2]]) {
  vdj_data_long <- vdj %>%
    gather(key = "gene_type", value = "gene", features)

  vdj_frequencies <- vdj_data_long %>%
    group_by(gene, .data[[ident_col]]) %>%
    summarise(freq = n()) %>%
    ungroup()

  count_matrix <- vdj_frequencies %>%
    spread(key = .data[[ident_col]], value = freq, fill = 0)
  count_matrix[is.na(count_matrix)] <- 0
  counts <- as.matrix(count_matrix %>% select(-gene))
  rownames(counts) <- count_matrix$gene
  coldata <- as.data.frame(ident_grp_table)
  coldata <- coldata[match(colnames(counts), coldata[[ident_col]]), ]
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                                design = as.formula(paste("~", grp_col)))
  dds <- dds[rowSums(counts(dds)) > 10, ]
  dds <- DESeq(dds, fitType="local")
  res <- results(dds)
  res$gene <- rownames(res)
  res$status <- ifelse(res$padj < 0.05 & res$log2FoldChange > 0,
                       paste0("Enriched in ", grp2),
                       ifelse(res$padj < 0.05 & res$log2FoldChange < 0,
                              paste0("Depleted in ", grp2), "Not Significant"))
  res
}

create_isotype_distribution_plot <- function(vgm,
                                             sel_celltype = NA,
                                             group_col = "clinical_condition",
                                             celltype_col = "celltype",
                                             isotype_cols = NA) {
  clinical_condition <- vgm[["VDJ"]][[group_col]]
  vdj_cgene <- vgm[["VDJ"]]$VDJ_cgene
  celltypes <- vgm[["VDJ"]][[celltype_col]]

  data <- data.frame(
    clinical_condition = clinical_condition,
    vdj_cgene = vdj_cgene,
    celltype = celltypes
  )

  if (!is.na(sel_celltype)) {
    data <- data %>% filter(celltype == sel_celltype)
  }

  data <- data %>%
    mutate(isotype = case_when(
      vdj_cgene == "IGHG1" ~ "IgG1",
      vdj_cgene == "IGHG2" ~ "IgG2",
      vdj_cgene == "IGHG3" ~ "IgG3",
      vdj_cgene == "IGHG4" ~ "IgG4",
      vdj_cgene == "IGHA1" ~ "IgA1",
      vdj_cgene == "IGHA2" ~ "IgA2",
      vdj_cgene == "IGHE" ~ "IgE",
      vdj_cgene == "IGHD" ~ "IgD",
      vdj_cgene == "IGHM" ~ "IgM",
      TRUE ~ "undet."
    ))

  data <- data %>% filter(isotype != "undet.")

  data_summary <- data %>%
    group_by(.data[[group_col]], isotype) %>%
    summarize(count = n(), .groups = "drop")

  chi_square_test <- chisq.test(table(data[[group_col]], data$isotype))

  significance_level <- ""
  if (chi_square_test$p.value < 0.001) {
    significance_level <- "***"
  } else if (chi_square_test$p.value < 0.01) {
    significance_level <- "**"
  } else if (chi_square_test$p.value < 0.05) {
    significance_level <- "*"
  } else {
    significance_level <- "-"
  }

  if (is.na(isotype_cols)) {
    isotype_cols <- c(
      "IgA1" = "#3d58a7",
      "IgA2" = "#f23325",
      "IgD" = "#76c161",
      "IgE" = "#e78ac3",
      "IgG1" = "#9553a1",
      "IgG2" = "#f7911e",
      "IgG3" = "#010202",
      "IgG4" = "#8da0cb",
      "IgM" = "#942e28",
      "undet." = "#7f7f7f")
  }

  plot_pie <- function(df, condition) {
    ggplot(df %>% filter(.data[[group_col]] == condition),
           aes(x = "", y = count, fill = isotype)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta = "y") +
      labs(title = condition, x = NULL, y = NULL) +
      theme_void() +
      scale_fill_manual(values = isotype_cols) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            plot.margin = unit(c(0, 0, 0, 0), "cm"),
            panel.spacing = unit(0, "lines"))
  }

  plots <- list()
  conditions <- unique(data[[group_col]])
  for (condition in conditions) {
    plots[[condition]] <- plot_pie(data_summary, condition)
  }

  legend_plot <- ggplot(data_summary, aes(x = "", y = count, fill = isotype)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    theme_void() +
    scale_fill_manual(values = isotype_cols) +
    theme(legend.position = "right")
  legend <- cowplot::get_legend(legend_plot)

  combined_pie_charts <- arrangeGrob(
    grobs = plots,
    ncol = length(conditions),
    padding = unit(0, "cm")
  )

  empty_plot <- ggplot() +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))

  significance_plot <- ggplot() +
    annotate("segment", x = 0.95, xend = 1.05, y = 1, yend = 1, size = 0.3) +
    annotate("text", x = 1, y = 0.875, label = significance_level, size = 6) +
    annotate("text", x = 1, y = 0.8, label = paste0("(χ² p-value: ",
                    formatC(chi_square_test$p.value,
                            format = "e", digits = 2), ")"), size = 2) +
    theme_void() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(panel.spacing = unit(0, "lines")) +
    ylim(0.5, 1.1)

  significance_combined <- arrangeGrob(
    empty_plot,
    significance_plot,
    empty_plot,
    ncol = 3,
    widths = unit(c(2.85,3,3.15), "null"),
    padding = unit(0, "cm")
  )

  combined_plots <- arrangeGrob(
    combined_pie_charts,
    empty_plot,
    significance_combined,
    ncol = 1,
    heights = unit(c(8, -0.5, 2), "cm"),
    padding = unit(0, "cm")
  )

  if (is.na(sel_celltype)) {
    title <- "Heavy Chain Isotype Distribution"
  } else {
    title <- paste0(sel_celltype, " Heavy Chain Isotype Distr.")
  }

  final_plot <- grid.arrange(
    combined_plots,
    legend,
    ncol = 2,
    widths = c(2.8, 0.5),
    top = textGrob(title, gp = gpar(fontsize = 16, fontface = "bold"))
  )

  final_plot
}

light_chain_Cgene_distribution_plot <- function(vgm, sel_celltype = NA) {
  # Extract data
  clinical_condition <- vgm[["VDJ"]]$clinical_condition
  vj_cgene <- vgm[["VDJ"]]$VJ_cgene
  celltypes <- vgm[["VDJ"]]$celltype

  data <- data.frame(
    clinical_condition = clinical_condition,
    vj_cgene = vj_cgene,
    celltype = celltypes
  )

  if (!is.na(sel_celltype)) {
    data <- data %>% filter(celltype == sel_celltype)
  }

  data <- data %>%
    mutate(chain = case_when(
      vj_cgene == "IGKC" ~ "Kappa",
      vj_cgene %in% c("IGLC1", "IGLC2", "IGLC3", "IGLC7") ~ "Lambda",
      TRUE ~ "Other"
    ))

  # if Kappa or Lambda does not exist, because it had value zero, create it
  # and add a small value to it, then filter out "Other" category
  if (!"Kappa" %in% data$chain) {
    data <- rbind(data, data.frame(clinical_condition = NA, vj_cgene = NA,
                                   celltype = NA, chain = "Kappa"))
  }
  if (!"Lambda" %in% data$chain) {
    data <- rbind(data, data.frame(clinical_condition = NA, vj_cgene = NA,
                                   celltype = NA, chain = "Lambda"))
  }
  data <- data %>% filter(chain %in% c("Kappa", "Lambda"))

  # Summarize the data for plotting
  data_summary <- data %>%
    group_by(clinical_condition, chain) %>%
    summarize(count = n())

  # Perform chi-square test
  # add a pseudo-value to data$chain if there is only one value

  chi_square_test <- chisq.test(table(data$clinical_condition, data$chain))

  # Determine significance level
  significance_level <- ""
  if (chi_square_test$p.value < 0.001) {
    significance_level <- "***"
  } else if (chi_square_test$p.value < 0.01) {
    significance_level <- "**"
  } else if (chi_square_test$p.value < 0.05) {
    significance_level <- "*"
  } else {
    significance_level <- "-"
  }

  # Define a function to create pie charts for each condition without legend
  plot_pie <- function(df, condition) {
    ggplot(df %>% filter(clinical_condition == condition),
           aes(x = "", y = count, fill = chain)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta = "y") +
      scale_fill_manual(values = c("Kappa" = "#008000", "Lambda" = "#ff8844")) +
      labs(title = condition, x = NULL, y = NULL) +
      theme_void() +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            plot.margin = unit(c(0, 0, 0, 0), "cm"),
            panel.spacing = unit(0, "lines"))
  }

  # Create the pie charts
  plots <- list()
  conditions <- unique(data$clinical_condition)
  for (condition in conditions) {
    plots[[condition]] <- plot_pie(data_summary, condition)
  }

  # Extract the legend from a separate plot
  legend_plot <- ggplot(data_summary, aes(x = "", y = count, fill = chain)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    theme_void() +
    scale_fill_manual(values = c("Kappa" = "#008000", "Lambda" = "#ff8844")) +
    theme(legend.position = "right")
  legend <- cowplot::get_legend(legend_plot)

  # Combine the pie charts without the legend
  combined_pie_charts <- arrangeGrob(
    grobs = plots,
    ncol = length(conditions),
    padding = unit(0, "cm")
  )

  empty_plot <- ggplot() +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))

  significance_plot <- ggplot() +
    annotate("segment", x = 0.95, xend = 1.05, y = 1, yend = 1, size = 0.3) +
    annotate("text", x = 1, y = 0.875, label = significance_level, size = 6) +
    annotate("text", x = 1, y = 0.8, label = paste0("(χ² p-value: ",
                              formatC(chi_square_test$p.value,
                                      format = "e", digits = 2), ")"),
             size = 2) +
    theme_void() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(panel.spacing = unit(0, "lines")) +
    ylim(0.5, 1.1)

  significance_combined <- arrangeGrob(
    empty_plot,
    significance_plot,
    empty_plot,
    ncol = 3,
    widths = unit(c(2.85,3,3.15), "null"),
    padding = unit(0, "cm")
  )

  # Combine the plots with an empty plot in between
  combined_plots <- arrangeGrob(
    combined_pie_charts,
    empty_plot,
    significance_combined,
    ncol = 1,
    heights = unit(c(8, -0.5, 2), "cm"),
    padding = unit(0, "cm")
  )

  if (is.na(sel_celltype)) {
    title <- "Light Chain Distribution"
  } else {
    title <- paste0(sel_celltype, " Light Chain Distr.")
  }

  # Arrange the combined plots and legend in a final layout
  final_plot <- grid.arrange(
    combined_plots,
    legend,
    ncol = 2,
    widths = c(2.8, 0.5),
    top = textGrob(title, gp = gpar(fontsize = 16, fontface = "bold"))
  )
  final_plot
}

gene_usage_signif <- function(vdj, celltype = NULL, celltype_column = NULL,
                              group_column = "clinical_condition",
                              genesets = c("VH", "JH", "VH_JH_pairing",
                                           "VL", "JL", "VL_JL_pairing",
                                           "VH_VL_pairing", "JH_JL_pairing",
                                           "VH_JH_VL_JL_pairing"),
                              n_permutations = 1000) {
  if (!is.null(celltype)) {
    if (is.null(celltype_column)) {
      stop("celltype_column must be specified when celltype is not NULL.")
    }
    message("Processing: ", celltype)
    vdj <- vdj[vdj[[celltype_column]] == celltype, ]
  } else {
    message("Processing: Total")
  }

  analyze_table <- function(tbl) {
    totals <- c(sum(tbl[, 1]), sum(tbl[, 2]))
    p_values <- apply(tbl, 1, function(counts) {
      chisq.test(matrix(counts, nrow = 2),
                 simulate.p.value = TRUE,
                 B = n_permutations)$p.value
    })
    p_adj <- p.adjust(p_values, method = "BH")
    stars <- ifelse(p_adj < 0.001, "***",
                    ifelse(p_adj < 0.01, "**",
                           ifelse(p_adj < 0.05, "*", "-")))

    calculate_log2fc <- function(counts, totals) {
      prop1 <- counts[1] / totals[1]
      prop2 <- counts[2] / totals[2]
      pseudo <- 1e-9
      log2((prop1 + pseudo) / (prop2 + pseudo))
    }
    log2fc_values <- apply(tbl, 1, calculate_log2fc, totals = totals)

    result <- data.frame(
      count1 = tbl[, 1],
      count2 = tbl[, 2],
      prop1 = tbl[, 1] / totals[1],
      prop2 = tbl[, 2] / totals[2],
      log2fc = log2fc_values,
      p_value = p_values,
      p_adj = p_adj,
      stars = stars
    )

    result <- result[order(-result$log2fc), ]
    return(result)
  }

  results_list <- list()

  if ("VH" %in% genesets) {
    results_list[["VH"]] <-
      analyze_table(table(vdj$VDJ_vgene, vdj[[group_column]]))
  }
  if ("JH" %in% genesets) {
    results_list[["JH"]] <-
      analyze_table(table(vdj$VDJ_jgene, vdj[[group_column]]))
  }
  if ("VH_JH_pairing" %in% genesets) {
    results_list[["VH_JH_pairing"]] <-
      analyze_table(table(paste(vdj$VDJ_vgene, vdj$VDJ_jgene, sep = "+"),
                          vdj[[group_column]]))
  }
  if ("VL" %in% genesets) {
    results_list[["VL"]] <-
      analyze_table(table(vdj$VJ_vgene, vdj[[group_column]]))
  }
  if ("JL" %in% genesets) {
    results_list[["JL"]] <-
      analyze_table(table(vdj$VJ_jgene, vdj[[group_column]]))
  }
  if ("VL_JL_pairing" %in% genesets) {
    results_list[["VL_JL_pairing"]] <-
      analyze_table(table(paste(vdj$VJ_vgene, vdj$VJ_jgene, sep = "+"),
                          vdj[[group_column]]))
  }
  if ("VH_VL_pairing" %in% genesets) {
    results_list[["VH_VL_pairing"]] <-
      analyze_table(table(paste(vdj$VDJ_vgene, vdj$VJ_vgene, sep = "+"),
                          vdj[[group_column]]))
  }
  if ("JH_JL_pairing" %in% genesets) {
    results_list[["JH_JL_pairing"]] <-
      analyze_table(table(paste(vdj$VDJ_jgene, vdj$VJ_jgene, sep = "+"),
                          vdj[[group_column]]))
  }
  if ("VH_JH_VL_JL_pairing" %in% genesets) {
    results_list[["VH_JH_VL_JL_pairing"]] <-
      analyze_table(table(paste(vdj$VDJ_vgene, vdj$VDJ_jgene,
                                vdj$VJ_vgene, vdj$VJ_jgene, sep = "+"),
                          vdj[[group_column]]))
  }

  results_df <- do.call(rbind, lapply(names(results_list), function(set_name) {
    if (nrow(results_list[[set_name]]) > 0) {
      data.frame(
        celltype = ifelse(is.null(celltype), "Total", celltype),
        tested_set = set_name,
        genes = rownames(results_list[[set_name]]),
        results_list[[set_name]],
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  }))

  return(results_df)
}

compute_clonotype_size_bins <- function(vgm) {
  clonotype_size_bin <- sapply(vgm[["GEX"]]$clonotype_size, function(x) {
    if (x == 1) {
      "1"
    } else if (x <= 5) {
      "2-5"
    } else if (x <= 20) {
      "6-20"
    } else {
      ">20"
    }
  })

  vgm[["GEX"]]$clonotype_size_bin <-
    factor(clonotype_size_bin, levels = c("1", "2-5", "6-20", ">20"))
  vgm <- transfer_gex_to_vdj(vgm, "clonotype_size_bin")
  vgm
}

# fix the clonotype output of the VDJ_clonotype_v3_w_enclone function
#
# for all duplicated barcodes
# - barcode_duplicated_group: all lines with same barcode
# for each barcode_duplicated_group:
# - if the clonotype_frequency is 1 for all lines in the barcode_duplicated_group:
#   then remove all lines, except the line with the lowest clonotype_id
#   (lowest numeric suffix of the clonotype_id column)
# - if there is one line with clonotype_frequency > 1 in the barcode_duplicated_group:
#   then remove all lines, except that one
# - if there are more than one lines with clonotype_frequency > 1 in the barcode_duplicated_group:
#   keep only the lowest clonotype id; also, for each other clonotype id in that barcode_duplicated_group,
#   change the clonotype_id of all lines outside of the barcode_duplicated_group with that clonotype_id
#   to the lowest clonotype_id
fix_enclone_clonotyping_output <- function(vdj_after) {
  duplicated_barcodes <- vdj_after$barcode[duplicated(vdj_after$barcode)]
  db_df <- subset(vdj_after, barcode %in% duplicated_barcodes)
  # split duplicated_barcodes by barcode:
  db_df_s <- split(db_df, db_df$barcode)
  clonotype_to_rm = c()
  # iterate over the db_df_s splits:
  for (barcode_duplicated_group in db_df_s) {
    # print barcode
    if (all(barcode_duplicated_group$clonotype_frequency == 1)) {
      lowest_clonotype_id <- paste0("clonotype", min(as.numeric(sub("clonotype", "", barcode_duplicated_group$clonotype_id))))
      clonotype_to_rm <- c(clonotype_to_rm,
                           barcode_duplicated_group$clonotype_id[barcode_duplicated_group$clonotype_id != lowest_clonotype_id])
    } else {
      # if there is a single line with clonotype_frequency > 1 in the barcode_duplicated_group:
      if (sum(barcode_duplicated_group$clonotype_frequency > 1) == 1) {
        clonotype_id_to_keep = barcode_duplicated_group$clonotype_id[barcode_duplicated_group$clonotype_frequency > 1]
        clonotype_to_rm <- c(clonotype_to_rm,
                             barcode_duplicated_group$clonotype_id[barcode_duplicated_group$clonotype_id != clonotype_id_to_keep])
      } else {
        # if there are more than one lines with clonotype_frequency > 1 in the barcode_duplicated_group:
        # keep only the lowest clonotype id; also, for each other clonotype id in that barcode_duplicated_group,
        # change the clonotype_id of all lines outside of the barcode_duplicated_group with that clonotype_id
        # to the lowest clonotype_id
      lowest_clonotype_id <- paste0("clonotype", min(as.numeric(sub("clonotype", "", barcode_duplicated_group$clonotype_id))))
      clonotype_to_rm <- c(clonotype_to_rm,
                           barcode_duplicated_group$clonotype_id[barcode_duplicated_group$clonotype_id != lowest_clonotype_id])
        # for each clonotype_id, identify the barcodes of all lines OUTSIDE of the barcode_duplicated_group
        # with that clonotype_id
        for (other_clonotype_id in barcode_duplicated_group$clonotype_id) {
          if (other_clonotype_id != lowest_clonotype_id) {
            print(paste0("changing: ", other_clonotype_id))
            barcodes_to_assign_to_clonotype_id <- vdj_after$barcode[vdj_after$clonotype_id == other_clonotype_id]
            barcodes_to_assign_to_clonotype_id <-
              barcodes_to_assign_to_clonotype_id[!barcodes_to_assign_to_clonotype_id %in% barcode_duplicated_group$barcode]
            vdj_after$clonotype_id[vdj_after$barcode %in% barcodes_to_assign_to_clonotype_id] <- lowest_clonotype_id
          }
        }
        print(clonotype_to_rm)
      }
    }
  }
  vdj_after <- vdj_after[!vdj_after$clonotype_id %in% clonotype_to_rm,]
}

