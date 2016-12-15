# Author: Walter Xie
# Accessed on 10 Mar 2016

#' Install all dependent packages.
#' 
#' @export
#' @examples 
#' installAllPackages()
installAllPackages <- function() {
  install.packages(c("vegan","vegetarian"))
  install.packages(c("grid", "gridExtra", "data.table", "xtable", "tools", "reshape2", "Matrix"))
  install.packages(c("gplots", "ggplot2", "RColorBrewer", "colorspace", "ggrepel", "NeatMap", "cowplot"))
  install.packages(c("ggdendro", "cluster", "clusterSim"))
  ###CommunityPhyloStru
  install.packages(c("ape", "picante", "plyr"))
  ###RDA
  install.packages(c("VIF", "scales", "corrgram", "pryr"))
  ###lasso
  install.packages(c("glmnet"))
}
