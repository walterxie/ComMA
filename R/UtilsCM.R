# Author: Walter Xie, Andrew Dopheide
# Accessed on 19 Nov 2015

#' Remove rows or colums from community matrix
#' 
#' Remove rows or colums in the matrix, whose sum of abundance is 
#' less and equal to the minimum abundance threshold. 
#' 
#' @param communityMatrix Community Matrix (OTU table), where rows are 
#' OTUs or individual species and columns are sites or samples. 
#' @param minAbund The minimum abundance threshold to remove rows/columns 
#' by row/column sum of abundance. For exampe, if minAbund=1, then remove 
#' all singletons appeared in only one sample. If minAbund=0, 
#' then remove all empty rows/columns. Default to 1 (singletons).
#' @param MARGIN 1 indicates rows, 2 indicates columns. Default to 1.
#' @param verbose More details. Default to TRUE.
#' @return the procceded community matrix
#' @keywords community matrix
#' @export
#' @examples 
#' remove singletons 
#' rmMinAbundance(communityMatrix, minAbund=1)
rmMinAbundance <- function(communityMatrix, minAbund=1, MARGIN=1, verbose=TRUE) {
  if (is.element(1, MARGIN)) {
    rm <- which(rowSums(communityMatrix)<=minAbund)
    if (length(rm)>0) 
      communityMatrix <- communityMatrix[-rm,]
    msg <- "from rows"
  } else if (is.element(2, MARGIN)) {
    rm <- which(colSums(communityMatrix)<=minAbund)
    if (length(rm)>0) 
      communityMatrix <- communityMatrix[,-rm]
    msg <- "from columns"
  }
  
  if(verbose) 
    cat("Remove", length(rm), msg, "having minimum abundance <=", minAbund, "!\n")
  
  communityMatrix
}

#' Transposed community matrix for \pkg{vegan} package
#' 
#' @param communityMatrix Community Matrix (OTU table)
#' @return the rotated community matrix
#' @keywords community matrix
#' @export
#' @examples communityMatrixT <- transposeCM(communityMatrix)
transposeCM <- function(communityMatrix) {
  if (!all(sapply(communityMatrix, is.numeric))) 
    stop("All community matrix elements have to be numeric type") 
  
  communityMatrixT <- as.data.frame(t(as.matrix(communityMatrix))) # transpose  
}



### Combine columns by sample name
mergeCMSample <- function(communityMatrix, sep) {
  if(missing(sep)) sep="-"
  
  colnames(communityMatrix) <- sapply(strsplit(colnames(communityMatrix), sep), "[[", 1) # Strip subplot letter by sep
  communityMatrix1 <- data.frame(matrix(ncol = 0, nrow = nrow(communityMatrix))) # Empty data.frame with required number of rows
  for(col in unique(colnames(communityMatrix))){
    cols <- communityMatrix[grep(col, colnames(communityMatrix))] # Find each pair of subplot columns
    cols1 <- as.data.frame(rowSums(cols)) # Add subplot pair together
    colnames(cols1) <- col
    communityMatrix1 <- cbind(communityMatrix1, cols1)
  }
  
  return(communityMatrix1)
}

# rowThr, colThr: remove row or/and column sum <= rowThr, colThr in communityMatrix,
# if rowThr, colThr = 0, then remove empty rows or/and columns
# mostAbundThr: keep the most abundant (mostAbundThr) OTUs only 
# return preprocessed communityMatrix
preprocessCM <- function(communityMatrix, keepSingleton, rowThr, colThr, mostAbundThr) { 
  # keepSingleton, abundPercThr, transverse, verbose
  # args <- list(...)
  if(missing(keepSingleton)) keepSingleton=TRUE
  if(missing(rowThr)) rowThr=0
  if(missing(colThr)) colThr=0
  if(missing(mostAbundThr)) mostAbundThr=0
  
  cat("keepSingleton = ", keepSingleton, "; rowThr = ", rowThr, "; colThr = ", colThr, 
      "; mostAbundThr = ", mostAbundThr, "\n") 
  cat("Original community matrix : samples = ", ncol(communityMatrix), ", OTUs/taxa = ", nrow(communityMatrix), ".\n") 
  
  # singletons
  if (!keepSingleton) {
    singletons <- which(rowSums(communityMatrix)==1)
    communityMatrix <- communityMatrix[-singletons,]
    cat("Remove", length(singletons) ,"singletons !\n")
    rm(singletons)		
  }	
  
  # this must be in front of filter column/row to avoid empty column
  if (mostAbundThr > 0) {
    communityMatrix <- keepMostAbundantRows(communityMatrix, mostAbundThr=mostAbundThr)
  }  
  
  communityMatrix <- rmMinAbundance(communityMatrix, minAbund=0, MARGIN=2)
  # filter column first to avoid empty rows after columns remvoed
  communityMatrix <- rmMinAbundance(communityMatrix, minAbund=0, MARGIN=1)
  
  # summary
  cat("Processed community matrix : samples = ", ncol(communityMatrix), ", OTUs/taxa = ", nrow(communityMatrix), ".\n") 
  
  return(communityMatrix)
}



