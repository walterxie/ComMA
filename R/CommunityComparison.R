# Author: Andrew Dopheide, Walter Xie
# Accessed on 25 Nov 2016


#' @name comparison
#' @title Compare community structure
#'
#' @description
#' Multivariate patterns of sample similarity were compared between pairwised community matrices 
#' using Procrustes and Mantel tests.
#' 
#' @details
#' \code{getDissimilarityList} calculates a list of similarity/dissimilarity 
#' given a list of community matrices and the distance metric.  
#' 
#' @param cm.list A list of community matrices.
#' @param metric,... The distance metric, default to "jaccard",
#' and other parametes passed to \code{\link{getDissimilarity}}.
#' @keywords comparison
#' @export
#' @examples 
#' jaccard.dist <- getDissimilarityList(cm.list, metric="jaccard")
#' jaccard.dist$dist.list
#' 
#' @rdname comparison
getDissimilarityList <- function(cm.list, metric="jaccard", ...){
  dist.list <- list()
  
  for (i in 1:length(cm.list)) {
    cm.name <- names(cm.list)[i]
    if (is.null(cm.name))
      cm.name <- paste0("community",i)
    cat("Calculate", metric, "dissimilarity for", cm.name, ".\n")
    
    dist.cm <- ComMA::getDissimilarity(cm.list[[i]], method=metric, ...)
    dist.list[[cm.name]] <- dist.cm
  }
  
  list(dist.list=dist.list, metric=metric)
}

#' @details
#' Mantel tests for comparisons between the overall community structure.  
#' 
#' @param dist.list A list of dissimilarity in \code{\link{dist}} 
#' between pairwised community matrices, 
#' which can be calculated by \code{\link{getDissimilarityList}}. 
#' @param method,permutations Parameters used in \code{\link{mantel}} in \pkg{vegan}.
#' @export
#' @examples 
#' corr.sign <- mantelComparison(jaccard.dist$dist.list)
#' 
#' @rdname comparison
mantelComparison <- function(dist.list, method="pearson", permutations = 999){
  corr.list <- list() # list for mantel results storage
  sign.list <- list()
  
  label.matched <- getLabelMatchedDist(dist.list)
  require(vegan)
  for(i in 1:length(label.matched$dist1)){
    man <- mantel(label.matched$dist1[[i]], label.matched$dist2[[i]], 
                  method=method, permutations = permutations)
    corr.list[[i]] <- man$statistic
    sign.list[[i]] <- man$signif
  }
  
  list(corr = corr.list, sign = sign.list, pairs=label.matched$pairs, 
       samples=label.matched$samples, same.samples=label.matched$same.samples)
}

#' @details
#' Procrustes comparisons between the overall community structure.  
#' 
#' @param scale,symmetric,permutations Parameters used in 
#' \code{\link{procrustes}} and \code{\link{protest}}.
#' @export
#' @examples 
#' procrustes <- procrustesComparison(jaccard.dist$dist.list)
#' 
#' @rdname comparison
procrustesComparison <- function(dist.list, scale = TRUE, symmetric = TRUE, permutations = 999){  
  proc.list <- list() # list to store procrustes objects (for plotting)
  prot.ss <- list() # list to store protest results 
  prot.corr <- list()
  prot.sign <- list()
  
  label.matched <- getLabelMatchedDist(dist.list)
  require(vegan)
  for(i in 1:length(label.matched$dist1)){
    # Generate MDS plots
    d1.mds <- metaMDS(label.matched$dist1[[i]])
    d2.mds <- metaMDS(label.matched$dist2[[i]])
    proc <- procrustes(d1.mds, d2.mds, scale = scale, symmetric = symmetric)
    prot <- protest(d1.mds, d2.mds, scale = scale, symmetric = symmetric, 
                    permutations = permutations)
    proc.list[[i]] <- proc
    prot.ss[[i]] <- prot$ss
    prot.corr[[i]] <- prot$t0
    prot.sign[[i]] <- prot$signif
  }

  list(proc = proc.list, prot.ss = prot.ss, prot.corr = prot.corr, prot.sign = prot.sign,
       pairs=label.matched$pairs, samples=label.matched$samples, same.samples=label.matched$same.samples)
}


###### Internal #####

getLabelMatchedDist <- function(dist.list) {
  dist1.list <- list()
  dist2.list <- list()
  pair.list <- list()
  samples.list <- list()
  same.samples.list <- list()
  x <- combn(dist.list, 2, simplify = F) # Get all pairs of dist matrices
  for(i in 1:length(x)){
    cat("Paired dist matrices : ", labels(x[[i]]), "\n")
    if( !all(labels(x[[i]][[1]]) == labels(x[[i]][[2]])) ) { # Different numbers of samples
      cat("Making samples match...\n")
      # Subset to shared samples
      keep <- intersect(labels(x[[i]][[1]]), labels(x[[i]][[2]]))
      m1 <- as.matrix(x[[i]][[1]])
      m2 <- as.matrix(x[[i]][[2]])
      d1 <- as.dist(m1[keep, keep])
      d2 <- as.dist(m2[keep, keep])
      same.samples.list[[i]] <- FALSE
    } else { # Same numbers of samples
      d1 <- x[[i]][[1]]
      d2 <- x[[i]][[2]]
      same.samples.list[[i]] <- TRUE
    }
    cat("Samples match ... ", all(names(d1) == names(d2)), "\n")
    dist1.list[[i]] <- d1
    dist2.list[[i]] <- d2
    pair.list[[i]] <- labels(x[[i]])
    samples.list[[i]] <- labels(d1)
  }
  list(dist1=dist1.list, dist2=dist2.list, pairs=pair.list, 
       samples=samples.list, same.samples=same.samples.list)
}



