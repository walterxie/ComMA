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
#' dist.list <- getDissimilarityList(cm.list, metric="jaccard")
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
#' @keywords classification
#' @export
#' @examples 
#' corr.sign <- mantelComparison(dist.list)
#' 
#' @rdname comparison
mantelComparison <- function(dist.list, method="pearson", permutations = 999){
  corr.sign <- list() # list for mantel results storage
  x <- combn(dist.list, 2, simplify = F) # Get all pairs of dist matrices
  for(i in 1:length(x)){
    cat("Paired dist matrices : ", labels(x[[i]]), "\n")
    if(length(labels(x[[i]][[1]])) != length(labels(x[[i]][[2]]))){ # Different numbers of samples
      cat("Making samples match...\n")
      # Subset to shared samples
      keep <- intersect(labels(x[[i]][[1]]), labels(x[[i]][[2]]))
      m1 <- as.matrix(x[[i]][[1]])
      m2 <- as.matrix(x[[i]][[2]])
      d1 <- as.dist(m1[keep, keep])
      d2 <- as.dist(m2[keep, keep])
    } else { # Same numbers of samples
      cat("Samples match...\n")
      d1 <- x[[i]][[1]]
      d2 <- x[[i]][[2]]
    }
    require(vegan)
    man <- mantel(d1, d2, method=method, permutations = permutations)
    t1 <- labels(x[[i]])[[1]] 
    t2 <- labels(x[[i]])[[2]]
    corr.sign[[i]] <- list(t1 = t1, t2 = t2, corr = man$statistic, sign = man$signif)
  }
  return(corr.sign)
}



### Procrustes comparisons ###
#do.procrustes <- function(dist.list, metric){
procrustesComparison <- function(dist.list){  
  p.results <- list() # list to store procrustes results 
  pro.list <- list() # list to store procrustes objects (for plotting)
  n <- 1 # counter for plot output
  x <- combn(dist.list, 2, simplify = F) # Get all pairs of dist matrices
  for(i in 1:length(x)){
    print(labels(x[[i]]))
    if(length(labels(x[[i]][[1]])) != length(labels(x[[i]][[2]]))){ # Different numbers of samples
      print("Making samples match...")
      # Subset to shared samples
      keep <- intersect(labels(x[[i]][[1]]), labels(x[[i]][[2]]))
      m1 <- as.matrix(x[[i]][[1]])
      m2 <- as.matrix(x[[i]][[2]])
      d1 <- as.dist(m1[keep, keep])
      d2 <- as.dist(m2[keep, keep])
    } else { # Same numbers of samples
      print("Samples match...")
      d1 <- x[[i]][[1]]
      d2 <- x[[i]][[2]]
    }
    # Generate MDS plots
    d1.mds <- metaMDS(d1)
    d2.mds <- metaMDS(d2)
    proc <- procrustes(d1.mds, d2.mds, scale = TRUE, symmetric = TRUE)
    prot <- protest(d1.mds, d2.mds, scale = TRUE, symmetric = TRUE, permutations = 999)
    p.results[[n]] <- list(t1 = labels(x[[i]])[[1]], t2 = labels(x[[i]])[[2]], 
                           pro_ss = prot$ss, pro_corr = prot$t0, pro_sig = prot$signif) 
    pro.list[[n]] <- proc
    names(pro.list)[[n]] <- paste(labels(x[[i]])[[1]], "vs.", labels(x[[i]])[[2]])
    n <- n + 1
  }
  return(list(p.results, pro.list))
}

