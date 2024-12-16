plotHighestExprs4Seurat <- function(SeuratObj, nfeat2plot = 50){
  # Required libraries
  require(ggplot2)
  require(reshape2)
  require(RColorBrewer)
  require(ggsci)
  # extracting count matrix
  mat_exp <- as.matrix(SeuratObj@assays$RNA@counts)
  # calculate sum of expression (sum of UMI)
  ave_exprs <- rowSums2(mat_exp)
  # ordering vector according to the level of detected genes
  ord_feat <- order(ave_exprs, decreasing = TRUE)
  selfeatures <- head(ord_feat, nfeat2plot)
  # Subtracting the expression matrix
  sub_mat <- as.data.frame(mat_exp[selfeatures, ])
  sub_mat$feature <- factor(rownames(sub_mat), levels = rev(rownames(sub_mat)))
  # converting dataframe into long format
  sub_mat_l <- melt(sub_mat, id.vars = "feature")
  # setting colors for the plot
  my_col <- colorRampPalette(rev(pal_futurama("planetexpress")(12)))(nfeat2plot)
  # make the plot
  my_plot <- ggplot(sub_mat_l, aes(x = value, y = feature))+
    geom_boxplot(aes(fill = feature), alpha = .75)+
    scale_fill_manual(values = my_col)+
    theme_light()+
    guides(fill=FALSE)+
    ylab("counts")
  my_plot
}
