# Author: Walter Xie, Andrew Dopheide
# Accessed on 20 Apr 2016

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

#' @details \code{transposeCM} returns a transposed community matrix for \pkg{vegan} package
#' 
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

#' @details \code{summaryCM.Vector} return a vector of summary of the 
#' community matrix, where \code{community.matrix} can be one column only.
#' The vector is c("reads","OTUs","Shannon","samples","singletons","doubletons").
#' 
#' @param digits The digits to \code{\link{round}} decimal places 
#' if number is not interger. Default to 2.
#' @keywords community matrix
#' @export
#' @examples 
#' summary.cm.vector <- summaryCM.Vector(community.matrix)
#'
#' @rdname utilsCM
summaryCM.Vector <- function(community.matrix, digits=2) {
  require(vegetarian)
  samples <- ncol(community.matrix)
  otus <- nrow(community.matrix)
  reads <- sum(community.matrix)
  rs <- rowSums(community.matrix)
  singletons <- sum(rs==1)
  doubletons <- sum(rs==2)
  shannon <- d(community.matrix,lev="gamma",q=1)
  
  return(prettyNum(c(reads, otus, round(shannon, digits), 
           samples, singletons, doubletons)))
}

#' @details \code{summaryCM} summarizes the community matrix.
#' 
#' @param has.total If 0, then only return abudence by samples (columns) of community matrix. 
#' If 1, then only return toal abudence. If 2, then return abudence by samples (columns) and total. 
#' Default to 1.
#' @param most.abund The threshold to define the number of the most abundent OTUs.
#' @keywords community matrix
#' @export
#' @examples 
#' summary.cm <- summaryCM(community.matrix)
#'
#' @rdname utilsCM
summaryCM <- function(community.matrix, most.abund, has.total=1, digits=2) {
  summary.cm <- data.frame(row.names = c("reads","OTUs","Shannon","samples","singletons","doubletons"))
  if (has.total!=1) {
    for (col.name in colnames(community.matrix)) 
      summary.cm[,col.name] <- summaryCM.Vector(community.matrix[,col.name], digits=digits)
  }
  if (has.total > 0) {
    summary.cm[,"total"] <- summaryCM.Vector(community.matrix, digits=digits)
    
    if (!missing(most.abund)) {
      if (most.abund > nrow(community.matrix))
        most.abund <- nrow(community.matrix)
      cat("Set most abundent OTUs threshold =", most.abund, ".\n")
      
      community.matrix <- community.matrix[order(rs, decreasing=TRUE),]
      cm <- community.matrix[1:most.abund,]
      col.name <- paste0("most.abund.", most.abund, ".otus)")
      
      summary.cm[,col.name] <- summaryCM.Vector(cm, digits=digits)
    }
  }
  return(summary.cm)
}


subsetCM <- function(cm.taxa, unclassified=0) {
  
}

#' @details \code{mostAbundantRows} trims data frame 
#' to the rows having most abundance given a threshold. 
#' The trimmed data frame having most abundant rows, 
#' such as community matrix of 150 most abundant OTUs.
#' 
#' @param most.abund The threshold to define the number 
#' of the most abundent OTUs. Default to 150.
#' @keywords community matrix
#' @export
#' @examples 
#' community.matrix <- getCommunityMatrix("16S", isPlot=TRUE, minAbund=1)
#' OTU100 <- mostAbundantRows(community.matrix, mostAbundant=100)
#'
#' @rdname utilsCM
mostAbundantRows <- function(community.matrix, most.abund=150) {
  if (nrow(community.matrix) < most.abund) 
    most.abund <- nrow(community.matrix)
  
  cat("Trim matrix to", most.abund, "rows having most abundance.\n") 
  
  rs <- rowSums(community.matrix)
  # order row sums decreasing
  ord<-order(rs, decreasing=TRUE) 
  community.matrix <- community.matrix[ord,]    
  
  return(community.matrix[1:most.abund,])
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
