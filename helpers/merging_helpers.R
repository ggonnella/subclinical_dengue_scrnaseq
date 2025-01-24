#
# Note: code of this function originally written by S.Mella
#
prepare_so_merging <- function(so, verbose=F){
  if (verbose) cat("Changing default assay to 'RNA'\n")
  DefaultAssay(so) <- "RNA"
  if (verbose) cat("Removing SCT assay\n")
  so@assays$SCT <- NULL
  so@meta.data <-
    so@meta.data[-grep(pattern = "SCT", colnames(so@meta.data))]
  if (verbose) cat("Removing graphs\n")
  for(ng in names(so@graphs)){
    so[[ng]] <- NULL
  }
  for(r in Seurat::Reductions(so)){
    so[[r]] <- NULL
  }
  if (verbose) cat("\t",unique(so$orig.ident, " done.\n"))
  return(so)
}
