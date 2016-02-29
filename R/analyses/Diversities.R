
# Author: Walter Xie, Alexei Drummond
# Accessed on 9 Sep 2015

library(vegetarian)

######## Jost diversity section #######
# communityMatrixT = t(communityMatrix), row is sample

#' Community data from file as a matrix where rows are OTUs or individual species and columns are sites or samples. 
#' Matrix elements are abundance data (e.g. counts, percent cover estimates).
#' @param communityMatrixT is abundances augment in d{vegetarian}, which is a transposed community matrix from file
diversity.df <- function(communityMatrixT=communityMatrixT) { 
  diversity.df <- data.frame(row.names=c("gamma", "alpha", "beta"))
  
  diversity.df$'q=0' <- c(
    d(communityMatrixT,lev="gamma",q=0),
    d(communityMatrixT,lev="alpha",q=0),
    d(communityMatrixT,lev="beta",q=0))
  
  diversity.df$'q=1' <- c(
    d(communityMatrixT,lev="gamma",q=1),
    d(communityMatrixT,lev="alpha",q=1),
    d(communityMatrixT,lev="beta",q=1))
  
  diversity.df$'q=2' <- c(
    d(communityMatrixT,lev="gamma",q=2),
    d(communityMatrixT,lev="alpha",q=2),
    d(communityMatrixT,lev="beta", q=2))
  
  colnames(diversity.df) <- c("$q=0$", "$q=1$", "$q=2$")
  rownames(diversity.df) <- c("$D_\\gamma(q)$", "$D_\\alpha(q)$", "$D_\\beta(q)$")
  
  return(diversity.df)
}

#' abundance (reads, gamme0) per sample
#' return 1-column data frame
abundancePerSample <- function(communityMatrixT, hasTotal) {
  if(missing(hasTotal)) hasTotal=TRUE
  
  # gamme0
  perSample <- data.frame(abundance=rowSums(communityMatrixT), stringsAsFactors=FALSE)
  rownames(perSample) <- rownames(communityMatrixT)
  
  if (hasTotal) {
    perSample <- rbind(perSample, sum(perSample[,1]))
    rownames(perSample)[nrow(perSample)] <- "Total"
  }
  
  return(perSample)
}

#' richness (OTUs/species) per sample
#' return 1-column data frame
richnessPerSample <- function(communityMatrixT, hasTotal) {
  if(missing(hasTotal)) hasTotal=TRUE
  
  # richness
  perSample <- data.frame(richness=rowSums(communityMatrixT > 0), stringsAsFactors=FALSE)
  rownames(perSample) <- rownames(communityMatrixT)
  
  if (hasTotal) {
    perSample <- rbind(perSample, sum(perSample[,1]))
    rownames(perSample)[nrow(perSample)] <- "Total"
  }
  
  return(perSample)
}

#' Shannon index (gamma1) per sample
#' return 1-column data frame
shannonPerSample <- function(communityMatrixT, digits = 2) {
  # Shannon
  #  gamma1 <- function(r) d(r,lev="gamma",q=1)
  #  perSample <- data.frame(Shannon=apply(communityMatrixT, 1, gamma1), stringsAsFactors=FALSE)
  perSample <- data.frame(row.names=rownames(communityMatrixT), stringsAsFactors=FALSE)
  for (i in 1:nrow(communityMatrixT)) 
    perSample[i,1] <- d(communityMatrixT[i,],lev="gamma",q=1)
  
  colnames(perSample)[1] <- "Shannon"
  
  return(round(perSample, digits))
}

######## alpha1 #######
#' effective alpha per sample
#' return one column matrix
alpha1 <- function(communityMatrixT) {    
  # including diagonal
  m.alpha1 <- matrix(0,nrow=nrow(communityMatrixT),ncol=1)	
  rownames(m.alpha1) <- c(rownames(communityMatrixT))
  for(i in 1:nrow(communityMatrixT)){				
    m.alpha1[i,1] <- d(communityMatrixT[i,],lev="gamma",q=1)				
  }
  
  return (m.alpha1) # 1 col matrix
}

######## Pair-wise turnovers #######

#' return a matrix cols and rows are sample names 
#' @param communityMatrixT a transposed community matrix for \pkg{vegetarian} 
#' @param diss.fun similarity/dissimilarity index, values are jaccard, horn.morisita, bray.curtis, and beta1-1. 
#' @param printProgressBar TRUE/FALSE, if missing, it will be nrow(row.pairs)>100
calculateDissimilarityMatrix <- function(communityMatrixT, diss.fun="beta1-1", printProgressBar) {    
  # including diagonal
  diss.matrix <- matrix(0,nrow=nrow(communityMatrixT),ncol=nrow(communityMatrixT))
  colnames(diss.matrix) <- c(rownames(communityMatrixT))
  rownames(diss.matrix) <- c(rownames(communityMatrixT))
  # row.pairs : each row is a pair of row number of communityMatrixT
  row.pairs <- t(combn(nrow(communityMatrixT),2))
  
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
      diss.matrix[row.pairs[n,2], row.pairs[n,1]] <- vegdist(communityMatrixT[row.pairs[n,],], method="jaccard", binary=TRUE)
    } else if (diss.fun=="horn.morisita") {
      # Horn-Morisita
      diss.matrix[row.pairs[n,2], row.pairs[n,1]] <- vegdist(communityMatrixT[row.pairs[n,],], method="horn", binary=FALSE)
    } else if (diss.fun=="bray.curtis") {
      # Bray-Curtis
      diss.matrix[row.pairs[n,2], row.pairs[n,1]] <- vegdist(communityMatrixT[row.pairs[n,],])
    } else { # diss.fun="beta1-1"
      # beta1-1
      diss.matrix[row.pairs[n,2], row.pairs[n,1]] <- d(communityMatrixT[row.pairs[n,],],lev="beta",q=1)-1
    }
  }
  if (printProgressBar) close(pb)
  
  return (diss.matrix)
}

#' Returns a distance matrix composed of pair-wise turnovers
#' Depends on: vegetarian library 
#' @ communityMatrixT a transposed community matrix for \pkg{vegetarian} 
TurnoverDist<-function(communityMatrixT){   
  to.table<-matrix(0,nrow=nrow(communityMatrixT),ncol=nrow(communityMatrixT))
  for(i in 1:nrow(communityMatrixT)){
    for(j in 1:nrow(communityMatrixT)){
      # For numerous communities of equal weights, the numbers equivalent of 
      # the Shannon beta diversity and the number of samples (N) can be used to 
      # calculate the turnover rate per sample (Equation 25 from Jost 2007, Harrison et al. 1992)
      to.table[i,j]<-turnover(communityMatrixT[c(i,j),])
    }
  }
  
  d <- as.dist(to.table)
  attr(d, "Labels") <- dimnames(communityMatrixT)[[1L]]
  
  return(d)
}


# COMPUTE TURNOVER TABLE
#
#d.turnover <- TurnoverDist(communityMatrixT)
#

# COMPUTER HORN-MORISITA OVERLAPS

#library(vegan)
#d.hornMorisita <- vegdist(communityMatrixT, method="horn", binary=FALSE)
#d.brayBin <- vegdist(communityMatrixT, method="bray", binary=TRUE)


#library(untb)
#cm_counts <- count(colSums(communityMatrixT))
#theta <- round(optimal.theta(cm_counts),2)
