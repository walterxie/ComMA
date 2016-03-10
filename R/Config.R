# Author: Walter Xie
# Accessed on 10 Mar 2016

#' Install all depdent packages.
#' 
#' @export
#' @examples 
#' installAllPackages()
installAllPackages <- function() {
  install.packages(c("vegan","vegetarian"))
  install.packages(c("grid", "gridExtra", "data.table", "xtable", "tools", "reshape2", "Matrix"))
  install.packages(c("gplots", "ggplot2", "RColorBrewer", "colorspace"))
  install.packages(c("ggdendro", "cluster"))
  ###CommunityPhyloStru
  install.packages(c("ape", "picante", "plyr"))
  ###RDA
  install.packages(c("VIF", "scales"))
}