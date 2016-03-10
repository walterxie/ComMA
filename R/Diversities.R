# Author: Walter Xie, Alexei Drummond
# Accessed on 9 Sep 2015

######## Jost diversity section #######
# t.communityMatrix = t(communityMatrix), row is sample

#' @name JostDiversity
#' @title Jost diversity from \pkg{vegetarian} Package
#'
#' @description Data input \strong{t.communityMatrix} is 
#' a transposed matrix of community matrix we defined in \pkg{ComMA}.
#' Community matrix from file is a matrix where rows are OTUs or individual species 
#' and columns are sites or samples. 
#' Matrix elements are abundance data (e.g. counts, percent cover estimates).
#' 
#' @param t.communityMatrix is abundances argument in \pkg{vegetarian} \code{\link{d}}, 
#' which is a transposed matrix of community matrix, 
#' where rows are plots (Use plots instead of subplots.), columns are OTUs.
#' @return 
#' \code{diversityTable} returns a 3x3 data frame: columns are levels of diversity c("gamma", "alpha", "beta"), 
#' rows are orders of the diversity measure c(0, 1, 2). For example,
#' \tabular{rrrr}{
#'    \tab $q=0$ \tab $q=1$ \tab $q=1$\cr
#'   $D_\\gamma(q)$ \tab 13922.000000 \tab 2501.693162 \tab 601.509610\cr
#'   $D_\\alpha(q)$ \tab 2238.392857 \tab 880.944977 \tab 251.127187\cr
#'   $D_\\beta(q)$ \tab 6.219641 \tab 2.839784 \tab 2.395239 
#' }
#' @export
#' @keywords diversity
#' @examples 
#' diversity.table <- diversityTable(t.communityMatrix)
#' 
#' @rdname JostDiversity
diversityTable <- function(t.communityMatrix) { 
  require(vegetarian)
  diversity.df <- data.frame(row.names=c("gamma", "alpha", "beta"))
  
  diversity.df$'q=0' <- c(
    d(t.communityMatrix,lev="gamma",q=0),
    d(t.communityMatrix,lev="alpha",q=0),
    d(t.communityMatrix,lev="beta",q=0))
  
  diversity.df$'q=1' <- c(
    d(t.communityMatrix,lev="gamma",q=1),
    d(t.communityMatrix,lev="alpha",q=1),
    d(t.communityMatrix,lev="beta",q=1))
  
  diversity.df$'q=2' <- c(
    d(t.communityMatrix,lev="gamma",q=2),
    d(t.communityMatrix,lev="alpha",q=2),
    d(t.communityMatrix,lev="beta", q=2))
  
  colnames(diversity.df) <- c("$q=0$", "$q=1$", "$q=2$")
  rownames(diversity.df) <- c("$D_\\gamma(q)$", "$D_\\alpha(q)$", "$D_\\beta(q)$")
  
  return(diversity.df)
}

#' abundance (reads, gamme0) per sample
#' return 1-column data frame
abundancePerSample <- function(t.communityMatrix, hasTotal=TRUE) {
  # gamme0
  perSample <- data.frame(abundance=rowSums(t.communityMatrix), stringsAsFactors=FALSE)
  rownames(perSample) <- rownames(t.communityMatrix)
  
  if (hasTotal) {
    perSample <- rbind(perSample, sum(perSample[,1]))
    rownames(perSample)[nrow(perSample)] <- "Total"
  }
  
  return(perSample)
}

#' richness (OTUs/species) per sample
#' return 1-column data frame
richnessPerSample <- function(t.communityMatrix, hasTotal=TRUE) {
  # richness
  perSample <- data.frame(richness=rowSums(t.communityMatrix > 0), stringsAsFactors=FALSE)
  rownames(perSample) <- rownames(t.communityMatrix)
  
  if (hasTotal) {
    perSample <- rbind(perSample, sum(perSample[,1]))
    rownames(perSample)[nrow(perSample)] <- "Total"
  }
  
  return(perSample)
}

#' Shannon index (gamma1) per sample
#' return 1-column data frame
shannonPerSample <- function(t.communityMatrix, digits = 2) {
  # Shannon
  #  gamma1 <- function(r) d(r,lev="gamma",q=1)
  #  perSample <- data.frame(Shannon=apply(t.communityMatrix, 1, gamma1), stringsAsFactors=FALSE)
  perSample <- data.frame(row.names=rownames(t.communityMatrix), stringsAsFactors=FALSE)
  for (i in 1:nrow(t.communityMatrix)) 
    perSample[i,1] <- d(t.communityMatrix[i,],lev="gamma",q=1)
  
  colnames(perSample)[1] <- "Shannon"
  
  return(round(perSample, digits))
}

######## alpha1 #######
#' effective alpha per sample
#' return one column matrix
alpha1 <- function(t.communityMatrix) {    
  # including diagonal
  m.alpha1 <- matrix(0,nrow=nrow(t.communityMatrix),ncol=1)	
  rownames(m.alpha1) <- c(rownames(t.communityMatrix))
  for(i in 1:nrow(t.communityMatrix)){				
    m.alpha1[i,1] <- d(t.communityMatrix[i,],lev="gamma",q=1)				
  }
  
  return (m.alpha1) # 1 col matrix
}

######## Pair-wise turnovers #######

#' Calculate similarity/dissimilarity distance matrix between samples.
#' 
#' @param diss.fun Similarity/dissimilarity index, values are jaccard, horn.morisita, 
#' bray.curtis, and beta1-1. Default to beta1-1, but it is slower than other indices.
#' @param printProgressBar Whether to print progress bar, if missing, it will be nrow(row.pairs)>100
#' @return 
#' \code{calculateDissimilarityMatrix} returns a \code{\link{matrix}} of similarity/dissimilarity, 
#' whose rows and columns are sample names
#' @export
#' @keywords diversity
#' @examples 
#' diss.matrix <- calculateDissimilarityMatrix(t.communityMatrix, diss.fun="jaccard")
#' 
#' @rdname JostDiversity
calculateDissimilarityMatrix <- function(t.communityMatrix, diss.fun="beta1-1", printProgressBar) {    
  # including diagonal
  diss.matrix <- matrix(0,nrow=nrow(t.communityMatrix),ncol=nrow(t.communityMatrix))
  colnames(diss.matrix) <- c(rownames(t.communityMatrix))
  rownames(diss.matrix) <- c(rownames(t.communityMatrix))
  # row.pairs : each row is a pair of row number of t.communityMatrix
  row.pairs <- t(combn(nrow(t.communityMatrix),2))
  
  cat("\nCalculating", diss.fun, "from", nrow(row.pairs), "pairs of samples.\n")
  
  if (missing(printProgressBar)) printProgressBar=nrow(row.pairs)>100
  if (printProgressBar) {
    flush.console()
    pb <- txtProgressBar(min=1, max=nrow(row.pairs), style = 3)
  }
  
  for (n in 1:nrow(row.pairs)) {
    if (printProgressBar) setTxtProgressBar(pb, n)
    if (diss.fun=="jaccard") {
      # Jaccard
      diss.matrix[row.pairs[n,2], row.pairs[n,1]] <- vegdist(t.communityMatrix[row.pairs[n,],], method="jaccard", binary=TRUE)
    } else if (diss.fun=="horn.morisita") {
      # Horn-Morisita
      diss.matrix[row.pairs[n,2], row.pairs[n,1]] <- vegdist(t.communityMatrix[row.pairs[n,],], method="horn", binary=FALSE)
    } else if (diss.fun=="bray.curtis") {
      # Bray-Curtis
      diss.matrix[row.pairs[n,2], row.pairs[n,1]] <- vegdist(t.communityMatrix[row.pairs[n,],])
    } else { # diss.fun="beta1-1"
      # beta1-1
      diss.matrix[row.pairs[n,2], row.pairs[n,1]] <- d(t.communityMatrix[row.pairs[n,],],lev="beta",q=1)-1
    }
  }
  if (printProgressBar) close(pb)
  
  return (diss.matrix)
}

#' Calculate pair-wise turnovers between samples.
#' 
#' @return 
#' \code{TurnoverDist} returns a \code{\link{dist}} composed of pair-wise turnovers.
#' @export
#' @keywords diversity
#' @examples 
#' turnover.dist <- TurnoverDist(t.communityMatrix)
#' 
#' @rdname JostDiversity
TurnoverDist<-function(t.communityMatrix){ 
  require(vegetarian)
  to.table<-matrix(0,nrow=nrow(t.communityMatrix),ncol=nrow(t.communityMatrix))
  for(i in 1:nrow(t.communityMatrix)){
    for(j in 1:nrow(t.communityMatrix)){
      # For numerous communities of equal weights, the numbers equivalent of 
      # the Shannon beta diversity and the number of samples (N) can be used to 
      # calculate the turnover rate per sample (Equation 25 from Jost 2007, Harrison et al. 1992)
      to.table[i,j]<-turnover(t.communityMatrix[c(i,j),])
    }
  }
  
  d <- as.dist(to.table)
  attr(d, "Labels") <- dimnames(t.communityMatrix)[[1L]]
  
  return(d)
}


# COMPUTE TURNOVER TABLE
#
#d.turnover <- TurnoverDist(t.communityMatrix)
#

# COMPUTER HORN-MORISITA OVERLAPS

#library(vegan)
#d.hornMorisita <- vegdist(t.communityMatrix, method="horn", binary=FALSE)
#d.brayBin <- vegdist(t.communityMatrix, method="bray", binary=TRUE)


#library(untb)
#cm_counts <- count(colSums(t.communityMatrix))
#theta <- round(optimal.theta(cm_counts),2)