# Author: Walter Xie
# Accessed on 1 June 2016

#' Use \code{\link{readHTMLTable}} to covert the report of code coverage 
#' into a list of two data frames or matrices of percentages.
#' 
#' @param doc HTML document which can be a file name or a URL or 
#' an already parsed HTMLInternalDocument, 
#' or an HTML node of class XMLInternalElementNode, 
#' or a character vector containing the HTML content to parse and process. 
#' @param ... Other arguments passed to \code{\link{readHTMLTable}}.
#' @return 
#' A list of two data frames or matrices, the 1st is Package, 2nd is Class. 
#' @keywords IO
#' @export
#' @examples 
#' code.coverage <- readCodeCoverage("code_coverage/index.html")
#' # $Package
#' #     Package Class % Method % Line %
#' # 1    all classes    38.2     26.1   26.6
#' # $Class
#' #     Package Class % Method % Line %
#' # 1    beast.app    36.4      8.9      1
#' # 2    beast.app.beastapp      20      2.7    0.5
readCodeCoverage <- function(doc, ...) {
  suppressMessages(suppressWarnings(require(XML)))
  tables <- readHTMLTable(doc)
  if (length(tables) != 2)
    stop("Invaild html of code coverage, tables", length(tables), "!= 2 !")
  names(tables) <- c("Package", "Class")
  
  for (ta in 1:length(tables)) {
    if (ncol(tables[[ta]]) != 4)
      stop("Invaild html of code coverage, table", tables[[ta]][1,1], ", column", ncol(tables[[ta]]), "!= 4 !")
    
    colnames(tables[[ta]]) <- c("Package", "Class %", "Method %", "Line %")
    for (col in 2:4) 
      tables[[ta]][,col] <- gsub("([0-9])%.*", "\\1", tables[[ta]][,col])
  }
  
  return(tables)
}
