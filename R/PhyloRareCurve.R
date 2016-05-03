


#' Rarefaction curve of Phylogenetic Diversity
#'
#' Use \code{\link{phylocurve}} in David Nipperess's package \pkg{PDcalc} 
#' to build rarefaction curves of phylogenetic diversity of multiple samples. 
#' The detail to see \url{https://github.com/davidnipperess/PDcalc} or 
#' \url{http://davidnipperess.blogspot.co.nz}.
#' 
#' Data input \strong{t.community.matrix} is 
#' a transposed matrix from community matrix we defined in \pkg{ComMA}.
#' \strong{phylo.tree} is a rooted tree of phylo object, 
#' which can get from \pkg{ape} \code{\link{read.tree}}. 
#' 
#' @param t.community.matrix A transposed matrix from community matrix, 
#' where rows are samples, columns are OTUs.
#' @param phylo.tree A rooted tree of phylo object.
#' @param min.size return NULL, if \code{min(rowSums(t.community.matrix)) < min.size}. Default to 100.
#' @param sample.sizes A vector of sample sizes to build phylo rare curve. 
#' If NULL, then create automatically from 1 to \code{max(rowSums(t.community.matrix))}.
#' @param points.1.min,points.min.max If \code{sample.sizes} is NULL, the number of data points is defined 
#' between 1 and \code{min(rowSums(t.community.matrix))} or and \code{max(rowSums(t.community.matrix))}. 
#' @param log.fn The CSV file name to store phylo rare data frame. If NULL, then do not log the file.
#' @param verbose Default to TRUE
#' @export
#' @examples 
#' phylo.rare.df <- getPhylorareDF(t.community.matrix, phylo.tree)
#' phylo.rare.df
getPhylorareDF <- function(t.community.matrix, phylo.tree, min.size=100, sample.sizes=NULL, 
                           points.1.min=6, points.min.max=9, log.fn="phylo-rare.csv", verbose=T) {
  rare.max <- max(rowSums(t.community.matrix))
  rare.min <- min(rowSums(t.community.matrix))
  
  if (rare.min < min.size) {
    warning(paste("Min sample size in the community", rare.min, "< threshold ", min.size, ", return NULL !\n"))
    return(NULL)
  }
  
  if (is.null(sample.sizes))
    sample.sizes <- c(round(exp(seq(log(1), log(rare.min), length.out = points.1.min)), digits = 0)[1:(points.1.min-1)], 
                      round(exp(seq(log(rare.min), log(rare.max), length.out = points.min.max)), digits = 0))
  if (verbose)
    cat("Phylo rare: subsampling by individual.\n") 
  
  require(PDcalc)
  for (ss in sample.sizes) {
    if (verbose)
      cat("Subsampling size =", ss, ".\n") 
    
    # individual (default), site or species
    phylo.rare <- PDcalc::phylocurve(t.community.matrix, phylo.tree, ss, subsampling = "individual", replace =F)
    # create phylo.rare.df
    if (which(sample.sizes == ss) == 1) 
      phylo.rare.df <- data.frame(row.names=rownames(phylo.rare), check.names=FALSE)
    
    if (! all(tolower(rownames(phylo.rare.df)) == tolower(rownames(phylo.rare))) )
      stop("Sample names do not match between phylo.rare and community matrix !")
    
    phylo.rare.df[,paste0("size.", ss)] <- phylo.rare[,1]
  }
  
  if (! is.null(log.fn)) {
    if (verbose)
      cat("Log phylo.rare.df to", log.fn, ".\n") 
    
    write.csv(phylo.rare.df, log.fn, row.names=TRUE, quote=FALSE)
  }
  
  return(phylo.rare.df)
}

