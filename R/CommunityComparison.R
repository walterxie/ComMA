# Author: Andrew Dopheide, Walter Xie
# Accessed on 25 Nov 2016


#' @name CommunityComparison
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
#' @rdname CommunityComparison
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
#' \code{mantelComparison} makes Mantel tests for comparisons 
#' between the overall community structure.  
#' 
#' @param dist.list A list of dissimilarity in \code{\link{dist}} 
#' between pairwised community matrices, 
#' which can be calculated by \code{\link{getDissimilarityList}}. 
#' @param method,permutations Parameters used in \code{\link{mantel}} in \pkg{vegan}.
#' @export
#' @examples 
#' mantel <- mantelComparison(jaccard.dist$dist.list)
#' 
#' @rdname CommunityComparison
mantelComparison <- function(dist.list, method="pearson", permutations = 999){
  m.list <- list() # list for mantel results storage

  label.matched <- getLabelMatchedDist(dist.list)
  require(vegan)
  for(i in 1:length(label.matched$dist1)){
    man <- mantel(label.matched$dist1[[i]], label.matched$dist2[[i]], 
                  method=method, permutations = permutations)
    l1 = label.matched$pairs[[i]][[1]] 
    l2 = label.matched$pairs[[i]][[2]]
    m.list[[i]] <- list(l1 = l1, l2 = l2, corr = man$statistic, sign = man$signif)
    cat("Pair (", l1, ",", l2, ") correlation =", man$statistic, ", significance =", man$signif, "\n")
  }
  m.result <- do.call("rbind", lapply(m.list, data.frame))
  
  list(m.df=m.result, m.list=m.list, pairs=label.matched$pairs, samples=label.matched$samples, 
       same.samples=label.matched$same.samples)
}

#' @details
#' \code{procrustesComparison} makes Procrustes comparisons 
#' between the overall community structure.  
#' 
#' @param scale,symmetric,permutations Parameters used in 
#' \code{\link{procrustes}} and \code{\link{protest}}.
#' @export
#' @examples 
#' procrustes <- procrustesComparison(jaccard.dist$dist.list)
#' 
#' @rdname CommunityComparison
procrustesComparison <- function(dist.list, scale = TRUE, symmetric = TRUE, permutations = 999){  
  proc.list <- list() # list to store procrustes objects (for plotting)
  prot.list <- list()
  
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
    l1 = label.matched$pairs[[i]][[1]]
    l2 = label.matched$pairs[[i]][[2]]
    prot.list[[i]] <- list(l1 = l1, l2 = l2, ss=prot$ss, corr=prot$t0, sign=prot$signif)
    cat("Pair (", l1, ",", l2, ") correlation =", prot$t0, ", significance =", prot$signif, 
        ", sum of squares ss =", prot$ss,"\n")
  }
  p.result <- do.call("rbind", lapply(prot.list, data.frame))
  
  list(proc = proc.list, prot.df=p.result, prot.list=prot.list, pairs=label.matched$pairs, 
       samples=label.matched$samples, same.samples=label.matched$same.samples)
}

#' @details
#' \code{getCorrTriMatrix} combines Mantel test (lower triangle) 
#' and Procrustes test (upper triangle) correlations into one matrix 
#' for comparisons.  
#' 
#' @param mantel.corr.list,prot.corr.list The list of correlations from  
#' Mantel test and Procrustes test.
#' @export
#' @examples 
#' corr.tri.m <- getCorrTriMatrix(mantel$m.df, procrustes$prot.df)
#' 
#' @rdname CommunityComparison
getCorrTriMatrix <- function(mantel.df, prot.df) {
  
  
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



