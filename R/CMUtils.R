# Author: Walter Xie, Andrew Dopheide
# Accessed on 20 Apr 2016
# 
# some tricks :
# aggregate(as.formula(paste(". ~", cols)), data=community.matrix, FUN=function(x) sum(x>0))
# df <- do.call(rbind, lapply(list, data.frame))

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
#' ComMA::rmMinAbundance(community.matrix, minAbund=2)
#'
#' @rdname utilsCM
rmMinAbundance <- function(community.matrix, minAbund=2, MARGIN=1, verbose=TRUE) {
  if (is.element(1, MARGIN)) {
    rm <- which(rowSums(community.matrix)<minAbund)
    if (length(rm)>0) 
      community.matrix <- community.matrix[-rm,]
    msg <- "rows"
  } else if (is.element(2, MARGIN)) {
    rm <- which(colSums(community.matrix)<minAbund)
    if (length(rm)>0) 
      community.matrix <- community.matrix[,-rm]
    msg <- "columns"
  }
  
  if(verbose) 
    cat("Remove", length(rm), msg, "whose minimum abundance <", minAbund, "!\n")
  
  community.matrix
}

#' @details \code{transposeDF} returns a transposed data frame, 
#' such as transposed community matrix for \pkg{vegan} package.
#' 
#' @keywords community matrix
#' @export
#' @examples 
#' t.community.matrix <- transposeDF(community.matrix)
#'
#' @rdname utilsCM
transposeDF <- function(community.matrix, to.numeric=TRUE) {
  if (!all(sapply(community.matrix, is.numeric))) {
    if (to.numeric)
      community.matrix <- ComMA::convertType(community.matrix)
    else
      stop("All community matrix elements have to be numeric type") 
  }
  # transpose, community.matrix must be numeric
  t.community.matrix <- as.data.frame(t(as.matrix(community.matrix)))  
}

#' @details \code{preprocessCM} exclude any samples with 
#' excessively low abundance.
#' 
#' @param cm A community matrix not transposed,
#' Columns are samples.
#' @param rm.samples Remove specified samples in a vector, 
#' it can be a keyword shared in sample names.
#' The vector will convert to a string separated by '|' to multi-samples. 
#' Default to empty vector to do nothing.
#' @param min.abund,mean.abund.thr Exclude any samples with excessively 
#' low abundance. Defaul \code{min.abund=5, mean.abund.thr=0.025}.
#' The final threshold takes the maximun value of 
#' \code{max(min.abund, mean(colSums(cm))*mean.abund.thr)}. 
#' @export
#' @examples 
#' cm <- preprocessCM(cm, rm.samples=c("CM30b51","CM30b58"))
#' 
#' @rdname utilsCM
preprocessCM <- function(cm, rm.samples=c(), min.abund=5, mean.abund.thr=0.025) {
  # remove specified samples, it can be keywords.
  if (length(rm.samples) > 0) {
    rm <- paste(rm.samples, collapse = "|")
    n.samples <- ncol(cm)
    cm <- cm[, !grepl(rm, colnames(cm), ignore.case = T)]
    cat("Drop", n.samples-ncol(cm), "samples containing : ", paste(rm.samples, collapse = ","), "\n") 
  }
  # exclude any samples with excessively low abundance
  print(summary(colSums(cm)))
  if (min.abund > 0 || mean.abund.thr > 0) {
    n.samples <- ncol(cm)
    max.thr <- max(min.abund, mean(colSums(cm))*mean.abund.thr)
    cm <- cm[, colSums(cm) > max.thr]
    cm <- cm[rowSums(cm) > 0, ] # Exclude any empty col 
    cat("Drop", n.samples-nrow(cm), "samples with low abundance <= ", max.thr, "\n") 
  }
  return(cm) 
}

#' @details \code{preprocessEnv} subsets the enviornmental variables 
#' and make log transform to soil chemistry variables.
#' 
#' @param env The enviornmental meta-data, where rows are samples 
#' and columns are enviornmental variables.
#' @param sel.env.var The vector of selected environmental variables, 
#' which can be colnames(env) or their indices. 
#' Defaul to an empty vector to choose all variables.
#' @param log.var,log.base It normally needs log transform to soil chemistry variables.
#' Use \code{\link{plotCorrelations}} to visualize variables and determine 
#' whether log transform should be applied. Default to no log transform. 
#' @export
#' @examples 
#' env <- preprocessEnv(env, sel.env.var=c(4,5,8,9,14:22), log.var=c(5:8,9:11))
#' 
#' @rdname utilsCM
preprocessEnv <- function(env, rm.samples=c(), sel.env.var=c(), log.var=c(), log.base=2) {
  # remove specified samples, it can be keywords.
  if (length(rm.samples) > 0) {
    rm <- paste(rm.samples, collapse = "|")
    n.samples <- nrow(env)
    env <- env[!grepl(rm, rownames(env), ignore.case = T), ]
    cat("Drop", n.samples-nrow(env), "samples containing : ", paste(rm.samples, collapse = ","), "\n") 
  }
  # select environmental variables
  if (length(sel.env.var) > 0) {
    env <- env[, sel.env.var]
    cat("Select", length(sel.env.var), "environmental variables : ", paste(colnames(env), collapse = ","), "\n") 
  }
  # select environmental variables
  if (length(log.var) > 0) {
    env[,log.var] <- log(env[,log.var], log.base)
    env[env == "-Inf"] <- 0 # Replace inf with zero
    cat("Log transform", length(log.var), "variables : ", paste(colnames(env[,log.var]), collapse = ","), 
        "at base", log.base, "\n") 
  }
  return(env)
}


#' @details \code{spilt.df} spilt a data frame into chunks of data frames 
#' having equal rows/columns.
#' 
#' @param spilt.to The number of sub-data-frame to spilt. It must >= 2.
#' Default to 2.
#' @keywords community matrix
#' @export
#' @examples 
#' cm.list <- spiltCM(community.matrix)
#'
#' @rdname utilsCM
spilt.df <- function(community.matrix, spilt.to=2, MARGIN=1, verbose=TRUE) {
  if (spilt.to < 2)
    stop("Invalid spilt.to", spilt.to, "!")
  
  if (MARGIN != 1) {
    size <- ceiling(ncol(community.matrix) / 2)
    msg <- "columns"
  } else { 
    size <- ceiling(nrow(community.matrix) / 2)
    msg <- "rows"
  }
  if (verbose)
    cat("Split a data frame into", spilt.to, ", roughly", size, msg, "each.")
  
  if (size < 2) 
    return(community.matrix)
  
  cm.list <- list()
  i = 1
  if (MARGIN != 1) {
    cm.cols.list <- split(1:ncol(community.matrix), ceiling(seq_along(1:ncol(community.matrix))/size))
    for (cm.cols in cm.cols.list) {
      cm.list[[i]] <- community.matrix[,cm.cols]
      i=i+1
    }
  } else {
    cm.rows.list <- split(1:nrow(community.matrix), ceiling(seq_along(1:nrow(community.matrix))/size))
    for (cm.rows in cm.rows.list) {
      cm.list[[i]] <- community.matrix[cm.rows,]
      i=i+1
    }
  }
  
  return(cm.list)
}


#' @details \code{mostAbundantRows} takes the given number of 
#' most abundant rows (OTUs) from original community matrix 
#' to form a new matrix. The new matrix will sort by both 
#' \code{rowSums} and \code{colSums} in decreasing by default.
#' 
#' @param most.abund The threshold to define the number 
#' of the most abundent OTUs. Default to 150.
#' @param row.decreasing,col.decreasing Should the sort 
#' decreasing order of \code{colnames} or \code{colSums}
#' be TRUE? Refer to \code{\link{order}}. If NULL, do nothing.
#' Default to TRUE.
#' @keywords community matrix
#' @export
#' @examples 
#' community.matrix <- getCommunityMatrix("16S", isPlot=TRUE, minAbund=1)
#' OTU100 <- mostAbundantRows(community.matrix, most.abund=100)
#'
#' @rdname utilsCM
mostAbundantRows <- function(community.matrix, most.abund=150, 
                             row.decreasing=TRUE, col.decreasing=TRUE) {
  if (nrow(community.matrix) < most.abund) 
    most.abund <- nrow(community.matrix)
  
  cat("Take", most.abund, "most abundant rows.\n") 
  
  rs <- rowSums(community.matrix)
  cs <- colSums(community.matrix)
  if (! is.null(row.decreasing)) {
    # order row sums decreasing
    community.matrix <- community.matrix[order(rs, decreasing=row.decreasing),]   
  } 
  if (! is.null(col.decreasing)) {
    community.matrix <- community.matrix[,order(cs, decreasing=col.decreasing)]   
  } 

  return(community.matrix[1:most.abund,])
}


#' @details \code{cmYAcrossX} aggregates a community matrix to another abundance matrix 
#' to show the number of OTUs (y-axis) simultaneously appeared at the number of samples (x-axis). 
#' The 'samples' is the number of samlpes listed in sequence,
#' the 'OTUs' is the number of OTUs simultaneously appeared only in that number of samlpes, 
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
  suppressMessages(suppressWarnings(require(Matrix)))
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

