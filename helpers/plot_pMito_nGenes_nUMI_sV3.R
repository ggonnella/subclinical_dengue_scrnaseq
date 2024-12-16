plot_pMito_nGenes_nUMI_sV3 <- function(seuratObj, title){
  require(ggplot2)
  require(gridExtra)
  colors <- c("green1", "blue1", "orange1","red1")
  alph = .6 # defining transparencies
  p.pMito <- ggplot(seuratObj@meta.data, aes(x = levels(factor(orig.ident)), y = percent.mito))+
    geom_jitter(aes(color = percent.mito.cl), alpha = alph)+
    geom_violin(trim = T, width = .5, alpha = .5)+
    scale_color_manual(values = colors)+
    xlab("")+
    guides(color = F)+
    theme_minimal()
  
  p.nGene <- ggplot(seuratObj@meta.data, aes(x = levels(factor(orig.ident)), y = nFeature_RNA))+
    geom_jitter(aes(color = percent.mito.cl), alpha = alph)+
    geom_violin(trim = T, width = .5, alpha = .5)+
    scale_color_manual(values = colors)+
    xlab("")+
    guides(color = F)+
    ggtitle(title)+
    theme_minimal()
  
  p.nUMI <- ggplot(seuratObj@meta.data, aes(x = levels(factor(orig.ident)), y = nCount_RNA))+
    geom_jitter(aes(color = percent.mito.cl), alpha = alph)+
    geom_violin(trim = T, width = .5, alpha = .5)+
    scale_color_manual(values = colors)+
    xlab("")+
    guides(color = F)+
    theme_minimal()
  
  grid.arrange(p.pMito, p.nGene, p.nUMI, ncol = 3, nrow= 1)
}