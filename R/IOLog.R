# Author: Walter Xie
# Accessed on 18 May 2016

#' Reads MCMC log files and return a data frame
#' whose column names are parameters and row names are the number
#' of states at each sample.
#'
#' @param file The file to read/write.
#' @param rm.na.col If TRUE, then remove all columns with all
#' missing values (NA). Default to TRUE.
#' @param ... Other arguments passed to \code{\link{readFile}}.
#' @keywords IO
#' @export
#' @examples
#' mcmc.log <- readMCMCLog("data/star.beast.log")
readMCMCLog <- function(file, rm.na.col=TRUE, ...) {
  suppressMessages(require(ComMA))
  mcmc.log <- ComMA::readFile(file, comment.char = "#", msg.file="MCMC log",
                              msg.col="parameters", msg.row="samples", ...)
  
  if (rm.na.col) {
    col.names <- ComMA::trimAll(names(mcmc.log))
    # exclude all columns with all NA, such as BEAST (< 2.4.1) log
    no.empty.col <- col.names[sapply( mcmc.log, function(x) !all(is.na(x)) )]
    if (length(col.names) != length(no.empty.col)) {
      mcmc.log <- mcmc.log[,no.empty.col]
      warning("Remove ", length(col.names) - length(no.empty.col), " empty column(s) !\n",
              "There are ", length(no.empty.col), " parameters for analysis !\n")
    }
  }
  
  attr(mcmc.log,"file") <- file
  return(mcmc.log)
}
