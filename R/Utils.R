# Utils
# Author: Walter Xie
# Accessed on 11 Mar 2016

#' Move rows to the last rows of data frame given a column value 
#' 
#' @param data A data frame.
#' @param col The column, which can be a string (column name) 
#' or number (index of the column).  
#' @param regex Regular expression for \code{\link{grep}} to indetify rows. 
#' Default to "unclassified". Use "^string$" to perform exact match.
#' @param ignore.case Refer to \code{\link{grep}}.
#' @return
#' The data frame of rows moved to the last.
#' @keywords utils
#' @export
#' @examples 
#' 
#' # move all rows having "unclassified" phyla to the last rows
#' data <- mvRowsToLast(data, "phylum")
#' # move all rows having exctly mathced "Bacteria" in column 1 to the last rows
#' data <- mvRowsToLast(data, 1, regex="^Bacteria$")
mvRowsToLast <- function(data, col, regex="unclassified", ignore.case = TRUE) {
  id.match <- grep(regex, data[,col], ignore.case = ignore.case)
  if (length(id.match) > 0)
    data <- data[c(setdiff(1:nrow(data), id.match),id.match),]
  return(data)
}


#' \code{prettyNumbers} provide pretty numbers with comma separator 
#' to a given data frame \code{df}.
#' 
#' @param digits, integer indicating the number of decimal places, default o 2.
#' @param drop.0.tail If TRUE, then remove all 0's decimal, such as ".00"
#' 
#' @export
#' @examples 
#' prettyNumbers(df)
prettyNumbers <- function(df, digits = 2, drop.0.tail=TRUE) {
  df <- ComMA::convertType(df)
  df <- round(df, digits)
  df <- format(df, big.mark=",", scientific=F)
  if (drop.0.tail && digits > 0) {
    pattern <- paste0( "\\.", paste0(rep("0",digits), collapse="") )
    df <- ComMA::gusbDF(pattern, "", df)
  }
  return(df)
}

#' Generate coordinates for 2 clusters. 
#' 
#' @source \url{http://stackoverflow.com/questions/2397097/how-can-a-data-ellipse-be-superimposed-on-a-ggplot2-scatterplot}.
#' 
#' @param n The number of points. Default to 100.
#' @param seed An integer seed for \code{\link{set.seed}}. Default to 101.
#' @keywords utils
#' @export
#' @examples 
#' df.clusters <- random2Clusters()
random2Clusters <- function(n=100, seed=101) {
  #bootstrap
  set.seed(seed)
  x <- rnorm(n, mean=2)
  y <- 1.5 + 0.4*x + rnorm(n)
  df <- data.frame(x=x, y=y, group="A")
  x <- rnorm(n, mean=2)
  y <- 1.5*x + 0.4 + rnorm(n)
  df <- rbind(df, data.frame(x=x, y=y, group="B"))
}

#' normalize given vector. 
#' 
#' @source \url{https://stat.ethz.ch/pipermail/r-help//2012-October/336676.html}.
#' 
#' @param vec A vector.
#' @keywords utils
#' @export
#' @examples 
#' normed <- as.data.frame(lapply(mtcars, normalize))
#' lapply(normed, range)
normalize <- function(vec) {
  (vec - min(vec, na.rm=TRUE))/(max(vec,na.rm=TRUE) - min(vec, na.rm=TRUE))
}


#' Convert data frame columns to different type, 
#' default to numeric type for columns. 
#' But it assign NA to actual string.
#' 
#' @source Modified from 
#' \url{http://stackoverflow.com/questions/2288485/how-to-convert-a-data-frame-column-to-numeric-type}.
#' 
#' @param df A data frame.
#' @param MARGIN,FUN Refer to \code{\link{apply}}.
#' @keywords utils
#' @export
#' @examples 
#' df <- convertType(df)
#' df <- convertType(df, FUN=as.character)
convertType <- function(df, FUN=as.numeric, stringsAsFactors=FALSE, check.names=FALSE) {
  df1 <- data.frame(suppressWarnings(lapply(df, FUN)), stringsAsFactors=stringsAsFactors, check.names=check.names)
  rownames(df1) <- rownames(df)
  return(df1)
}

#' Return the first \code{n} elements. 
#' If \code{n} is greater than length,
#' then return the whole vector or list.
#' 
#' @param vect A Vector or List.
#' @param n The first n elements.
#' @keywords utils
#' @export
#' @examples 
#' first.n <- firstN(1:10, 3)
#' firstN(1:10, 20)
firstN <- function(vect, n) {
  vect[1:ifelse(length(vect) < n, length(vect), n)]
}



