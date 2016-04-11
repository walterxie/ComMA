#' ComMA R documentation
#'
#' This package is a collection of tools for community matrix analyse, 
#' especially for eDNA data, but it also can be used for traditional species survey data.
#' The package not only wraps bioinformatic, ecological methods and figure plotting  
#' into one-line functions, but also create utility functions to procced data easily. 
#' All functions are available by typing \code{help(package="ComMA")} in R. 
#'
#' \strong{Community matrix} in \pkg{ComMA} is a community data from file as a matrix,
#' where rows are OTUs (therefore this matrix is called \emph{OTU table}) or individual species 
#' and columns are sites or samples. 
#' Matrix elements are abundance data or proportion (e.g. counts, percentage). For example,
#' \tabular{rrrr}{
#'   OTU_id \tab plot01 \tab plot02\tab ...\cr
#'   OTU_1 \tab 1 \tab 0 \tab ...\cr
#'   OTU_2 \tab 100 \tab 200 \tab ...\cr
#'   OTU_3 \tab 56 \tab 3 \tab ...
#' }
#' 
#' Variable \strong{t.community.matrix} is a transposed matrix from community matrix 
#' we defined here, where columns are OTUs or individual species and rows are sites or samples. 
#' It is also the abundances argument in \pkg{vegetarian} \code{\link{d}}.
#' 
#' The \strong{taxa.table} is a data frame to contain taxonomic classifications of OTUs.  
#' Rows are OTUs or individual species, which have to be the subset of (or matching) rows 
#' in the community matrix from \code{getCommunityMatrix}). 
#' Columns are taxonomy at the rank, such as c("superkingdom", "kingdom", "phylum", 
#' "class", "order", "family", "genus", "species"), or full taxonomic lineage but it is optional. 
#' If it is from RDP, then column "confidence" is required.
#' For example,
#' \tabular{rrrr}{
#'   OTU_id \tab kingdom \tab phylum\tab ...\cr
#'   OTU_1 \tab Bacteria \tab Firmicutes \tab ...\cr
#'   OTU_2 \tab Bacteria \tab Proteobacteria \tab ...\cr
#'   OTU_3 \tab Not assigned \tab Not assigned \tab ...
#' }
#'  
#' The enviornmental data \strong{env} is also called as the meta data of samples, 
#' where rows are sites or samples, and columns are the measurement, such as 
#' elevation, tempurature, soil chemistry, forest type, etc.
#' For example,
#' \tabular{rrrr}{
#'   Plot \tab Elevation(m) \tab pH\tab ...\cr
#'   plot01 \tab 50 \tab 5.89 \tab ...\cr
#'   plot02 \tab 90 \tab 5.12 \tab ...\cr
#'   plot03 \tab 160 \tab 5.46 \tab ...
#' }
#' 
#' The \strong{phyloTree} is a rooted tree of \code{\link{phylo}} object, 
#' which is created from \pkg{ape} \code{\link{read.tree}}.
#' 
#' @source \url{https://github.com/walterxie/ComMA}
#' 
#' Install all depdent packages: \code{\link{installAllPackages()}} 
#' 
#' 
#' 
"_PACKAGE"
#> [1] "_PACKAGE"
