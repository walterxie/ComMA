# Author: Walter Xie, Andrew Dopheide, Alexei Drummond
# Accessed on 15 Sep 2016


#' @name distance
#' @title Calculate various multivariate distance metrics
#'
#' @description Data input \strong{cm} is community matrix defined in \pkg{\link{ComMA}}.
#' Rows are OTUs or individual species, and columns are sites or samples.  
#' It is a transposed matrix of \code{\link{vegdist}} input.

#' @details Calculate various multivariate distance metrics among samples.
#' 
#' @param cm The community matrix.
#' @param tree The \pkg{\link{ape}} tree object required to caculate \code{\link{UniFrac}}.
#' @param is.transposed If TRUE, then \code{cm} is already the transposed matrix 
#' of \code{\link{vegdist}} input. Default to FASLE.
#' @param method The method to calculate similarity/dissimilarity distance. 
#' The options are jaccard, horn.morisita, 
#' bray.curtis, beta1-1, and unwt.unif, wt.unif. 
#' The last two are unweighted or weighted \code{\link{UniFrac}} 
#' from \pkg{\link{phyloseq}}. 
#' Default to beta1-1, but it is slow.
#' @param progressBar Whether to print progress bar 
#' only when method="beta1-1", 
#' if TRUE, it will be nrow(row.pairs)>100.
#' @return 
#' \code{getDissimilarity} returns a \code{\link{matrix}} of similarity/dissimilarity, 
#' whose rows and columns are sample names
#' @export
#' @keywords distance
#' @examples 
#' dist.jaccard <- getDissimilarity(cm, method="jaccard")
#' 
#' @rdname distance
getDissimilarity <- function(cm, tree=NA, is.transposed=FALSE, method="beta1-1", progressBar=FALSE) {  
  # if is.transposed=F as default, then transpose to vegdist input
  if (!is.transposed)
    cm <- ComMA::transposeDF(cm)
  
  require(vegan)
  if (method=="jaccard") {
    # Jaccard
    return(vegdist(cm, method="jaccard", binary=TRUE))
  } else if (method=="horn.morisita") {
    # Horn-Morisita
    return(vegdist(cm, method="horn"))
  } else if (method=="bray.curtis") {
    # Bray-Curtis
    return(vegdist(cm, method = "bray"))
  } else if (method=="beta1-1") { 
    # method="beta1-1"
    return (ComMA::beta1minus1(cm, progressBar=progressBar))
  } else if(grepl("unif", method)){ # Unifrac distances
    require(phyloseq)
    physeq <- phyloseq(otu_table(cm, taxa_are_rows = FALSE), phy_tree(tree))

    if(method == "unwt.unif") # unweighted
      return( phyloseq::UniFrac(physeq, weighted = FALSE, normalized = TRUE) )
    else if(method == "wt.unif") # weighted
      return( phyloseq::UniFrac(physeq, weighted = TRUE, normalized = TRUE) )
  } 
  warning("Cannot find method ", method, "to calculate similarity/dissimilarity distance !\n")
  return (NULL)
}

# beta1-1
beta1minus1 <- function(cm, progressBar=FALSE, return.dist=TRUE) {
  # including diagonal
  dist.matrix <- matrix(0,nrow=nrow(cm),ncol=nrow(cm))
  colnames(dist.matrix) <- c(rownames(cm))
  rownames(dist.matrix) <- c(rownames(cm))
  # row.pairs : each row is a pair of row number of cm
  row.pairs <- t(combn(nrow(cm),2))
  
  cat("\nCalculating beta1-1 from", nrow(row.pairs), "pairs of samples.\n")
  if (progressBar) {
    flush.console()
    pb <- txtProgressBar(min=1, max=nrow(row.pairs), style = 3)
  }
  
  require(vegetarian)
  for (n in 1:nrow(row.pairs)) {
    if (progressBar) 
      setTxtProgressBar(pb, n)
    # beta1-1
    dist.matrix[row.pairs[n,2], row.pairs[n,1]] <- d(cm[row.pairs[n,],],lev="beta",q=1)-1
  }
  if (progressBar) 
    close(pb)
  if (return.dist)
    return(as.dist(dist.matrix))
  return(dist.matrix)
}


######## Pair-wise turnovers #######

#' @details Calculate pair-wise turnovers between samples.
#' 
#' @return 
#' \code{TurnoverDist} returns a \code{\link{dist}} composed of pair-wise turnovers.
#' @export
#' @keywords distance
#' @examples 
#' turnover.dist <- TurnoverDist(t.community.matrix)
#' 
#' @rdname distance
TurnoverDist<-function(t.community.matrix){ 
  require(vegetarian)
  turnover.table<-matrix(0,nrow=nrow(t.community.matrix),ncol=nrow(t.community.matrix))
  for(i in 1:nrow(t.community.matrix)){
    for(j in 1:nrow(t.community.matrix)){
      # For numerous communities of equal weights, the numbers equivalent of 
      # the Shannon beta diversity and the number of samples (N) can be used to 
      # calculate the turnover rate per sample (Equation 25 from Jost 2007, Harrison et al. 1992)
      turnover.table[i,j]<-turnover(t.community.matrix[c(i,j),])
    }
  }
  
  d <- as.dist(turnover.table)
  attr(d, "Labels") <- dimnames(t.community.matrix)[[1L]]
  
  return(d)
}

######## alpha1 #######
#' effective alpha per sample
#' return one column matrix
#alpha1 <- function(t.community.matrix) {    
# including diagonal
#  m.alpha1 <- matrix(0,nrow=nrow(t.community.matrix),ncol=1)	
#  rownames(m.alpha1) <- c(rownames(t.community.matrix))
#  for(i in 1:nrow(t.community.matrix)){				
#    m.alpha1[i,1] <- d(t.community.matrix[i,],lev="gamma",q=1)				
#  }
#  
#  return (m.alpha1) # 1 col matrix
#}

#library(untb)
#cm_counts <- count(colSums(t.community.matrix))
#theta <- round(optimal.theta(cm_counts),2)
