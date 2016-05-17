# Author: Walter Xie
# Accessed on 18 May 2016


#' @name IOLog
#' @title IO functions to read BEAST log
#'
#' @description IO functions to read BEAST 1 \url{http://beast.bio.ed.ac.uk} 
#' or BEAST 2 \url{http://www.beast2.org} log files.
#' 
#' @details 
#' \code{readBEASTLog} reads BEAST log files and return a data frame 
#' whose column names are parameters and row names are samples.
#' 
#' @param file The file to read/write.
#' @keywords IO
#' @export
#' @examples 
#' beast.log <- readBEASTLog("data-raw/star.beast.log")
#' 
#' @rdname IOLog
readBEASTLog <- function(file, verbose=TRUE) { 
  beast.log <- ComMA::readFile(file, verbose=verbose, msg.file="BEAST log", 
                               msg.col="parameters", msg.row="samples")
  
  attr(beast.log,"file") <- file
  return(beast.log)
}

