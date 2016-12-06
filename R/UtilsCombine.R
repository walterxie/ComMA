# Utils
# Author: Walter Xie
# Accessed on 29 Nov 2016

#' @name UtilsCombine
#' @title Utils to combine data frames or matrices into the required data format
#' 
#' @details 
#' \code{getTriMatrix} converts pairwised comparison result 
#' into a symmetric triangular matrix.
#' The pairwised comparison result is stored in a data frame.
#' Two of their columns must have pair names in the same order, 
#' which determine the matrix's row names and column names.
#' The 1st data frame goes to the lower triangle of the matrix, 
#' and the 2nd is in upper triangle.
#' 
#' @param df The data frame containing pairwised comparison result, 
#' such as correlations.  
#' For example, a data frame can be 
#' \tabular{rrrr}{
#'   l1 \tab l2 \tab corr\tab sign\cr
#'   16S \tab 18S \tab 0.827 \tab 0.001\cr
#'   16S \tab ITS \tab 0.585 \tab 0.001\cr
#'   18S \tab ITS \tab 0.729 \tab 0.001
#' }
#' @param key Two column names containing the pair names,
#' which determine the output matrix's row names and column names.
#' It must make row and column names same.
#' @param value The column's values are going to fill in the output matrix.
#' @param na.to.0 Logical, replace all NA to 0, as default.
#' @param order.by A vector to order matrix rows and columns. 
#' It must be a subset or same as their names.
#' @keywords utils
#' @export
#' @examples 
#' corr.tri <- getTriMatrix(corr.df)
#' 
#' @rdname UtilsCombine 
getTriMatrix <- function(df, row.col=c("l1", "l2"), value="corr", na.to.0=TRUE, order.by=c()) {
  if (!all(c(row.col, value) %in% colnames(df)))
    stop("Invalid inputs: cannot find column names ", paste(c(row.col, value), collapse = ","), " !")
  
  df2 <- df
  df2[,row.col[1]] <- df[,row.col[2]]
  df2[,row.col[2]] <- df[,row.col[1]]
  df2 <- rbind(df2, df) # Duplicating/reversing data makes symmetric matrix
  
  require(reshape2)
  m <- dcast(df2, as.formula(paste(row.col, collapse ="~")), value.var = value)
  rownames(m) <- m[,row.col[1]]
  # rm column l1
  m <- m[ , -which(colnames(m) %in% row.col)]

  # order rows and columns
  if (length(order.by) > 1) {
    if ( !all( is.element(order.by, rownames(m)) ) )
      stop("order.by must be a subset or same as row/colnames !")
    m <- m[match(order.by, rownames(m)), ]
    m <- m[, match(order.by, colnames(m))]
  } else {
    m <- m[order(rownames(m)), ]
    m <- m[, order(colnames(m))]
  }
  
  if (na.to.0)
    m[is.na(m)] <- 0
  
  if (!all(rownames(m) == colnames(m)))
    stop("Invaild result: key 'row.col' must make row/colnames same ! \n", 
         "Rownames: ", paste(rownames(m), collapse = ","), ". \n",
         "Colnames: ", paste(rownames(m), collapse = ","))
  return(m)
}

#' @details \code{combineTriMatrix} combines two symmetric triangular matrices into one.
#' The 1st matrix goes to the lower triangle in the combined matrix,
#' and the 2nd to upper triangle.
#' 
#' @param tri.m1,tri.m2 Two symmetric triangular matrices, 
#' which can be the output from \code{getTriMatrix}.
#' @keywords utils
#' @export
#' @examples 
#' corr.sign.tri <- combineTriMatrix(corr.tri, sign.tri)
#' 
#' @rdname UtilsCombine 
combineTriMatrix <- function(tri.m1, tri.m2) {
  if (!all(rownames(tri.m1) == rownames(tri.m2)) || !all(colnames(tri.m1) == colnames(tri.m2)))
    stop("Invalid inputs: row or column names not same !")
  
  tri.m1[upper.tri(tri.m1)] <- NA # lower tri 
  tri.m2[lower.tri(tri.m2)] <- NA # upper tri 
  m <- tri.m1
  m[upper.tri(m)] <- tri.m2[upper.tri(tri.m2)] # Combine matrices
  return(m)
}


#' @details \code{combineTwoDF} combines two data frames or matrices with a same structure in one, 
#' put all values in 2nd data frame into brackets.
#' 
#' @param dfm A data frame or matrix.
#' @param dfm2 The 2nd data frame or matrix whose values are into brackets.
#' @param rm.zero Default to TRUE to remove all " (0)".
#' @param return.df,... Default to TRUE to return a data frame, otherwise a matrix.
#' @keywords utils
#' @export
#' @examples 
#' df <- combineTwoDF(df, df2, stringsAsFactors=FALSE)
#' 
#' @rdname UtilsCombine 
combineTwoDF <- function(dfm, dfm2, rm.zero=TRUE, return.df=TRUE, ...) {
  if (nrow(dfm) != nrow(dfm2) || ncol(dfm) != ncol(dfm2)) 
    stop("Two data frames must have a same structure !")
  
  dfm <- as.matrix(dfm)
  dfm2 <- as.matrix(dfm2)
  dfm.comb <- matrix( paste0(trimSpace(dfm), " (", trimSpace(dfm2), ")"), 
          nrow=nrow(dfm), dimnames=dimnames(dfm) )
  
  if (rm.zero) 
    dfm.comb <- gsub(" \\(0\\)", "", dfm.comb)
  
  if (return.df)
    dfm.comb <- data.frame(dfm.comb, check.names=FALSE, ...)
  return(dfm.comb)
}

#' @details \code{mergeBy} is an improved function of 
#' \code{\link{merge}} \code{by} column(s).
#' 
#' @param x,y data frames, or objects to be coerced to one.
#' @param warning.msg Logical; if TRUE as default, then print warning message 
#' if any row is dropped after merge.
#' @param by,... Specifications of the column(s) used for merging, 
#' and more arguments passed to \code{\link{merge}}.
#' @param rm.by Logical; if TRUE, then remove the column(s) specified
#' by \code{by} after merge. Default to FALSE. 
#' But if \code{by="row.names"}, then copy the column 'Row.names' created 
#' by \code{\link{merge}} into rownames before drop it.
#' @keywords utils
#' @export
#' @examples 
#' df <- mergeBy(df, df2)
#' 
#' @rdname UtilsCombine 
mergeBy <- function(x, y, by="row.names", rm.by=FALSE, warning.msg=TRUE, ...) {
  xy <- merge(x, y, by=by, ...)
  
  if ( warning.msg && (nrow(xy) != nrow(x) || nrow(xy) != nrow(y)) ) 
    warning("Lossing rows after merge !  nrow(xy) = ", nrow(xy), 
            ", nrow(x) = ", nrow(x), ", nrow(y) = ", nrow(y), " !")
  
  if (rm.by) {
    # mv 'by' such as Row.names to rownames, if it is only 1 col 
    if (tolower(by[1]) == "row.names") {
      rownames(xy) <- xy$Row.names
      by="Row.names"
    } 
    xy <- xy[,-which(colnames(xy) %in% by)]
    cat("Drop column(s)", paste(by, collapse = ","), "after merge.\n")
  }

  return(xy)
}


#' @details \code{mergeListOfDF} merges a list of data frames into one 
#' using \code{\link{mergeBy}}.
#' 
#' @param df.list a list of data frames, or objects to be coerced to one.
#' @param suffixes The vector of suffixes added to distinguish colmun names 
#' merged by different data frames. Its length must equal to \code{df.list}.   
#' @keywords utils
#' @export
#' @examples 
#' df <- mergeListOfDF(df.list, suffixes=names(df.list))
#' df <- mergeListOfDF(df.list, by="gene", suffixes=c("16S","18S","ITS"))
#' 
#' @rdname UtilsCombine 
mergeListOfDF <- function(df.list, by="row.names", rm.by=TRUE, suffixes=c(), ...) {
  df.whole <- as.data.frame(df.list[[1]])
  
  if (length(suffixes) > 0) {
    if (length(suffixes) != length(df.list))
      stop("Invalid inputs : length(suffixes) != length(df.list) !")
    suffixes <- ComMA::trimSpace(suffixes, ".") # replace space to .
    cm.names <- suffixes
    # add suffix to 1st df cols 
    colnames(df.whole) <- paste(colnames(df.whole), suffixes[1], sep = ".")
  } else if (is.null(names(df.list))) {
    cm.names <- 1:length(df.list)
  } else {
    cm.names <- names(df.list)
  }
  cat("Load data 1 : ", cm.names[1], "\n")
  
  if (length(df.list) > 1) { # multiple cm
    for (i in 2:length(df.list)) {
      df.tmp <- as.data.frame(df.list[[i]])
      cat("Load data", i, ": ", cm.names[i], "\n")
      if (length(suffixes) > 0)
        colnames(df.tmp) <- paste(colnames(df.tmp), suffixes[i], sep = ".")
      # mv Row.names col to rownames
      df.whole <- ComMA::mergeBy(df.whole, df.tmp, by=by, rm.by=rm.by)
    }
  }
  return(df.whole)
}


