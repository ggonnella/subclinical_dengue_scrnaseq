factor_resolution <- function(SeuratObj){
  tmp_df <- SeuratObj@meta.data
  ind <- grep("res.", colnames(tmp_df))
  print(head(tmp_df[ind]))
  for(i in ind){
    tmp_df[,i] <- factor(tmp_df[,i])
    print(levels(tmp_df[,i]))
    tmp_df[,i] <- factor(tmp_df[,i], levels = 0 : (length(levels(tmp_df[,i])) - 1) )
    print(levels(tmp_df[,i]))
    cat("\n\n")
  }
  
  SeuratObj@meta.data <- tmp_df
  
  return(SeuratObj)
  
}
