remove_info <- function(seuratObj){
  cat("Changing default assay to 'RNA'\n")
  DefaultAssay(seuratObj) <- "RNA"
  cat("Removing SCT assay\n")
  seuratObj@assays$SCT <- NULL
  seuratObj@meta.data <-
    seuratObj@meta.data[-grep(pattern = "SCT", colnames(seuratObj@meta.data))]
  cat("Removing graphs\n")
  for(ng in names(seuratObj@graphs)){
    seuratObj[[ng]] <- NULL
  }
  for(r in Seurat::Reductions(seuratObj)){
    seuratObj[[r]] <- NULL
  }
  cat("\t",unique(seuratObj$orig.ident, " done.\n"))
  return(seuratObj)
}
