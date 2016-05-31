# Author: Walter Xie
# Accessed on 1 June 2016



#' Get lines containing Java method name.
#' The regular expression is modified from 
#' \url{http://stackoverflow.com/questions/68633/regex-that-will-match-a-java-method-declaration}.
#' 
#' @param linn A character vector of lines read from the class file. 
#' Refer to \code{\link{readLines}}.
#' @return 
#' A character vector of lines containing Java method name. 
#' @keywords IO
#' @export
#' @examples 
#' method.lines <- getMethodLines(linn)
getMethodLines <- function(linn) {
  pattern <- paste0("((public|private|protected|static|final|native|synchronized|abstract|transient|void)+\\s)",
                  "+[\\$_\\w\\<\\>\\[\\]]*\\s+[\\$_\\w]+\\([^\\)]*\\)?\\s*\\{?[^\\}]*\\}?")
  linn[grepl(pattern, linn, perl=T)]
}
