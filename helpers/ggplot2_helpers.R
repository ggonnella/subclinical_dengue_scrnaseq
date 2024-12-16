library(ggplot2)
library(grid)
split_print_legend <- function(plt) {
  print(plt + theme(legend.position = "none"))
  tmp <- ggplot_gtable(ggplot_build(plt)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  grid.newpage()
  grid.draw(legend) 
}