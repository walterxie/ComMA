# Author: Walter Xie, Andrew Dopheide
# Accessed on 19 Nov 2015

#' @name utilsCM
#' @title Utils to preprocess community matrix
#'
#' @description Utils to preprocess community matrix, 
#' such as removing OTUs by different filters, 
#' and aggregating matrix by different ways.
#' 
#' @details 
#' \code{rmMinAbundance} returns the subset matrix of given community matrix, 
#' by removing rows or colums whose sum of abundance is less than the minimum abundance threshold. 
#' 
#' @param community.matrix Community matrix (OTU table), where rows are 
#' OTUs or individual species and columns are sites or samples. See \code{\link{ComMA}}.
#' @param minAbund The minimum abundance threshold to remove rows/columns 
#' by row/column sum of abundance. For exampe, if minAbund=2, then remove 
#' all singletons appeared in only one sample. If minAbund=1, 
#' then remove all empty rows/columns. Default to 2 (singletons).
#' @param MARGIN 1 indicates rows, 2 indicates columns. Default to 1.
#' @param verbose More details. Default to TRUE.
#' @keywords community matrix
#' @export
#' @examples 
#' # remove singletons 
#' rmMinAbundance(community.matrix, minAbund=2)
#'
#' @rdname utilsCM
rmMinAbundance <- function(community.matrix, minAbund=2, MARGIN=1, verbose=TRUE) {
  if (is.element(1, MARGIN)) {
    rm <- which(rowSums(community.matrix)<minAbund)
    if (length(rm)>0) 
      community.matrix <- community.matrix[-rm,]
    msg <- "from rows"
  } else if (is.element(2, MARGIN)) {
    rm <- which(colSums(community.matrix)<minAbund)
    if (length(rm)>0) 
      community.matrix <- community.matrix[,-rm]
    msg <- "from columns"
  }
  
  if(verbose) 
    cat("Remove", length(rm), msg, "having minimum abundance <", minAbund, "!\n")
  
  community.matrix
}

#' @details \code{transposeCM} provides a transposed community matrix for \pkg{vegan} package
#' 
#' @return the rotated community matrix
#' @keywords community matrix
#' @export
#' @examples 
#' communityMatrixT <- transposeCM(community.matrix)
#'
#' @rdname utilsCM
transposeCM <- function(community.matrix) {
  if (!all(sapply(community.matrix, is.numeric))) 
    stop("All community matrix elements have to be numeric type") 
  
  communityMatrixT <- as.data.frame(t(as.matrix(community.matrix))) # transpose  
}

#' @details \code{cmYAcrossX} aggregates a community matrix to a data frame 
#' to show the number of OTUs (y-axis) across the number of samples (x-axis). 
#' The 'samples' is the number of samlpes in sequence,
#' the 'OTUs' is the number of OTUs appeared only in that number of samlpes, 
#' and the 'reads' is the number of reads assigned to those OTUs.
#' 
#' @param terms The terms to be used in names of the data frame, 
#' which will be shown in the graph if using \code{\link{ggBarYAcrossX}}. 
#' Default to c("OTUs", "samples", "reads"). 
#' Please be careful of the order if any change.
#' @keywords community matrix
#' @export
#' @examples 
#' community.matrix <- getCommunityMatrix("16S", isPlot=TRUE, minAbund=1)
#' cm.aggre <- cmYAcrossX(community.matrix)
#' print(cm.aggre, row.names = FALSE)
#'
#' @rdname utilsCM
cmYAcrossX <- function(community.matrix, terms=c("OTUs", "samples", "reads")) {
  require(Matrix)
  row.count.sum <- data.frame(row.names = rownames(community.matrix))
  row.count.sum[, terms[2]] <- apply(community.matrix, MARGIN=1, function(x) sum(x>0))
  row.count.sum$reads <- apply(community.matrix, MARGIN=1, sum)
  
  cm.aggre.c <- aggregate(as.formula(paste(". ~", terms[2])), data=row.count.sum, function(x) sum(x>0))
  names(cm.aggre.c)[names(cm.aggre.c)==terms[3]] <- terms[1] # 1: OTUs, 3: reads
  cm.aggre.s <- aggregate(as.formula(paste(". ~", terms[2])), data=row.count.sum, FUN=sum)
  cm.aggre <- merge(cm.aggre.c, cm.aggre.s, by = terms[2]) # 2: samples
  
  return(cm.aggre)
}



#' @details \code{sumColumns} sums the columns by the \emph{n}th substring defined in column names.
#' 
#' @param sep The seperator to get the \emph{n}th substring from column names. Default to dash '-'.
#' @param nth The \emph{n}th substring. Default to 1 (first).
#' @keywords community matrix
#' @export
#' @examples 
#' # by subpl
#' community.matrix <- getCommunityMatrix("16S", isPlot=FALSE, minAbund=1)
#' colSums(community.matrix)
#' communityMatrix1 <- sumColumns(community.matrix)
#' colSums(communityMatrix1)
#'
#' @rdname utilsCM
sumColumns <- function(community.matrix, sep="-", nth=1) {
  colnames(community.matrix) <- sapply(strsplit(colnames(community.matrix), split=sep), "[[", nth) # Strip subplot letter by sep
  communityMatrix1 <- data.frame(matrix(ncol = 0, nrow = nrow(community.matrix))) # Empty data.frame with required number of rows
  for(col in unique(colnames(community.matrix))){
    cols <- community.matrix[grep(col, colnames(community.matrix))] # Find each pair of subplot columns
    cols1 <- as.data.frame(rowSums(cols)) # Add subplot pair together
    colnames(cols1) <- col
    communityMatrix1 <- cbind(communityMatrix1, cols1)
  }
  
  return(communityMatrix1)
}

#' @details \code{mergeRowsByColumnValue} \code{\link{aggregate}}s 
#' given community matrix by applying a specified function (e.g. sum) 
#' to the values of given columns.
#' 
#' @param ... The community matrix column names for \code{mergeRowsByColumnValue}.
#' @param FUN The function to compute the summary statistics 
#' which can be applied to all data subsets. Refer to \code{\link{aggregate}}.
#' @keywords community matrix
#' @export
#' @examples 
#'
#' @rdname utilsCM
mergeRowsByColumnValue <- function(community.matrix, ..., FUN=sum) {
  cols <- paste(list(...), collapse="+")
  
  aggregate(as.formula(paste(". ~", cols)), data=community.matrix, FUN=FUN)
}

#' @details \code{mostAbundantRows} trims data frame 
#' to the rows having most abundance given a threshold. 
#' The trimmed data frame having most abundant rows, 
#' such as community matrix of 150 most abundant OTUs.
#' 
#' @param mostAbundant The threshold to return rows having most abundance. 
#' If \code{nrow(community.matrix) < mostAbundant}, then use \code{nrow}. Default to 150.
#' @keywords community matrix
#' @export
#' @examples 
#' community.matrix <- getCommunityMatrix("16S", isPlot=TRUE, minAbund=1)
#' OTU100 <- mostAbundantRows(community.matrix, mostAbundant=100)
#'
#' @rdname utilsCM
mostAbundantRows <- function(community.matrix, mostAbundant=150) {
  if (nrow(community.matrix) < mostAbundant) 
    mostAbundant <- nrow(community.matrix)
  
  cat("Trim matrix to", mostAbundant, "rows having most abundance.\n") 
  
  rs <- rowSums(community.matrix)
  # order row sums decreasing
  ord<-order(rs, decreasing=TRUE) 
  community.matrix <- community.matrix[ord,]    

  return(community.matrix[1:mostAbundant,])
}

# rowThr, colThr: remove row or/and column sum <= rowThr, colThr in community.matrix,
# if rowThr, colThr = 0, then remove empty rows or/and columns
# mostAbundThr: keep the most abundant (mostAbundThr) OTUs only 
# return preprocessed community.matrix
preprocessCM <- function(community.matrix, keepSingleton, rowThr, colThr, mostAbundThr) { 
  # keepSingleton, abundPercThr, transverse, verbose
  # args <- list(...)
  if(missing(keepSingleton)) keepSingleton=TRUE
  if(missing(rowThr)) rowThr=0
  if(missing(colThr)) colThr=0
  if(missing(mostAbundThr)) mostAbundThr=0
  
  cat("keepSingleton = ", keepSingleton, "; rowThr = ", rowThr, "; colThr = ", colThr, 
      "; mostAbundThr = ", mostAbundThr, "\n") 
  cat("Original community matrix : samples = ", ncol(community.matrix), ", OTUs/taxa = ", nrow(community.matrix), ".\n") 
  
  # singletons
  if (!keepSingleton) {
    singletons <- which(rowSums(community.matrix)==1)
    community.matrix <- community.matrix[-singletons,]
    cat("Remove", length(singletons) ,"singletons !\n")
    rm(singletons)		
  }	
  
  # this must be in front of filter column/row to avoid empty column
  if (mostAbundThr > 0) {
    community.matrix <- keepMostAbundantRows(community.matrix, mostAbundThr=mostAbundThr)
  }  
  
  community.matrix <- rmMinAbundance(community.matrix, minAbund=1, MARGIN=2)
  # filter column first to avoid empty rows after columns remvoed
  community.matrix <- rmMinAbundance(community.matrix, minAbund=1, MARGIN=1)
  
  # summary
  cat("Processed community matrix : samples = ", ncol(community.matrix), ", OTUs/taxa = ", nrow(community.matrix), ".\n") 
  
  return(community.matrix)
}



