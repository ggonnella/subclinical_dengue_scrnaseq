calculate_percentMito_sV3 <- function(seuratObj, mitoGenes, plot = TRUE, SC_gal_tools_path = NULL){
  if(is.null(SC_gal_tools_path)){
    cat("You need to indicate the directory where the function plot_pMito_nGenes_nUMI_sV3 is.\n")
  }
  else{
    source(paste0(SC_gal_tools_path,"plot_pMito_nGenes_nUMI_sV3.R"))
    title = substitute(seuratObj)
    percent.mito <- Matrix::colSums(seuratObj@assays$RNA@counts[rownames(seuratObj@assays$RNA@counts) %in% mitoGenes,]) / 
      Matrix::colSums(seuratObj@assays$RNA@counts)
    seuratObj@meta.data$percent.mito <- percent.mito
    seuratObj@meta.data$percent.mito.cl <- cut(seuratObj@meta.data$percent.mito, 
                                               breaks = c(-1, .1, .15, .2, 1))
    if(plot){plot_pMito_nGenes_nUMI_sV3(seuratObj = seuratObj, title)}
    return(seuratObj)
  }
}
