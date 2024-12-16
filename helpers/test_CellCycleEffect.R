test_CellCycleEffect <- function(seuratObj, seuratObj_Name, s_genes, g2m_genes, ret_so=FALSE){
  
  cat("\n### Processing ", seuratObj_Name, "\n")
  
  cat("\n##### Number of S genes in ", seuratObj_Name, " : \n")
  s_genes_table <- as.data.frame(table(s_genes %in% rownames(seuratObj)))
  
  cat("\n", paste( apply(s_genes_table,1, paste, collapse = " : "), collapse  = "\n\n"))
  
  cat("\n\n##### Number of G2M genes in ", seuratObj_Name, " : \n")
  g2m_genes_table <- as.data.frame(table(g2m_genes %in% rownames(seuratObj)))
  
  cat("\n", paste( apply(g2m_genes_table,1, paste, collapse = " : "), collapse  = "\n\n"))
  
  cat("\n\n#### Normalizing Seurat Object by nUMI count\n")
  seuratObj <- NormalizeData(seuratObj, assay = "RNA")
  
  cat("\n#### Score cells for cell cycle\n")
  seuratObj <- CellCycleScoring(object = seuratObj,
                                g2m.features = g2m_genes[g2m_genes %in% rownames(seuratObj)],
                                s.features = s_genes[s_genes %in% rownames(seuratObj)])
  
  cat("\n##### Frequency of the different cell cycle phases\n")
  ccp_df <- as.data.frame.array(table(seuratObj$Phase))
  ccp_df$phase <- rownames(ccp_df)
  
  cat("\n", paste( apply(ccp_df[,2:1],1, paste, collapse = " : "), collapse  = "\n\n"))
  
  cat("\n\n#### Preparing data for PCA\n")
  cat("\nIdentify most variable genes\n")
  seuratObj <- FindVariableFeatures(seuratObj, 
                                    selection.method = "vst",
                                    nfeatures = 2000, 
                                    verbose = FALSE,
                                    assay = "RNA")
  
  cat("\n#### Scaling the counts\n")
  seuratObj <- ScaleData(seuratObj, assay = "RNA", verbose = FALSE)
  
  cat("\n#### Perform PCA\n")
  seuratObj <- RunPCA(seuratObj, 
                      assay = "RNA", 
                      reduction.name = "pcaRNA", 
                      reduction.key = "pcaRNA_",
                      verbose = FALSE)
  
  cat("\n#### plot the PCA colored by cell cycle phase\n")
  ccg_colors <- RColorBrewer::brewer.pal(length(levels(seuratObj$Phase)), name = "Set1")
  
  if(ret_so){  
    plot(
      DimPlot(seuratObj,
              reduction = "pcaRNA",
              group.by= "Phase",
              cols = ccg_colors)
    )
    plot(
      DimPlot(seuratObj,
              reduction = "pcaRNA",
              group.by= "Phase",
              split.by = "Phase",
              cols = ccg_colors)
    )
    return(seuratObj)
  }
  else{
    plot(
      DimPlot(seuratObj,
              reduction = "pcaRNA",
              group.by= "Phase",
              split.by = "Phase",
              cols = ccg_colors)
    )
  }
  
}