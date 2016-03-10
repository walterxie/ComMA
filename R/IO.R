# Author: Walter Xie
# Accessed on 3 Mar 2016

#' Create folder from \code{\link{file.path}} if not exist.
#' 
#' @param dir.path The folder path constructed by \code{\link{file.path}}.
#' @export
#' @examples 
#' figDir <- file.path(workingPath, "figures")
#' mkdir(figDir) 
mkdir <- function(dir.path) {
  require(tools)
  if (!file.exists(dir.path)) {    
    dir.create(dir.path)    
  }
  cat("\nConfig : setup", dir.path, "\n")
}

#' Read a file to return a data frame. 
#' 
#' If the file extension is \emph{csv}, 
#' then use \code{\link{read.csv}}, otherwise use \code{\link{read.table}}.
#' 
#' @param file The 1st row is column names, the 1st column is row names.
#' @param sep Only used for non \emph{csv} file. Default to tab "\\t". 
#' @return 
#' A data frame from the file, such as
#' \tabular{rrrr}{
#'   OTU_id \tab plot01 \tab plot02\tab ...\cr
#'   OTU_1 \tab 1 \tab 0 \tab ...\cr
#'   OTU_2 \tab 100 \tab 200 \tab ...\cr
#'   OTU_3 \tab 56 \tab 3 \tab ...
#' }
#' @export
#' @examples 
#' communityMatrix <- readFile("16S.txt")
readFile <- function(file, sep="\t") { 
  require(tools)
  # sep="\t" only work for non csv file
  if (tolower(file_ext(file))=="csv") {
    df <- read.csv(file, head=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
  } else {
    # be careful read.table bug   
    df <- read.table(file, sep=sep, header=T, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)  
  }	
  return(df)
}

#' Write the data fram (table) to a file.
#' 
#' @param df A data frame. 
#' @param file If the file extension is \emph{csv}, 
#' then use \code{\link{write.csv}}, 
#' otherwise use \code{\link{write.table}}.
#' @export
#' @examples 
#' writeTable(df, file.path)
writeTable <- function(df, file){
  require(tools)
  if (tolower(file_ext(file))=="csv") {
    write.csv(df, file, quote=FALSE)
  } else { # .tsv .txt
    #write.table bug: mistake to start col names from the 1st cell
    cat("",colnames(df),file=file,sep="\t")
    cat("\n",file=file, append=TRUE)
    write.table(df, file, sep ="\t", quote=FALSE, col.names=FALSE, append=TRUE)
  }
}

#' Print data fram (table) as Latex format to either file or console.
#' 
#' @param df A data frame.
#' @param caption Latex table caption.
#' @param label Latex table label.
#' @param file If NULL, then print the results to console, otherwise print them to the file. Default to NULL. 
#' @param align Refer to \code{\link{xtable}}.
#' @param digits Refer to \code{\link{xtable}}. 
#' @export
#' @examples 
#' tableFile <- file.path(workingPath, "report.tex")
#' printXTable(data.frame, caption = "Phylogenetic beta diversity", 
#'             label = "tab:pd:beta", file=tableFile)
printXTable <- function(df, caption, label, file=NULL, align = NULL, digits = NULL) {
  require(xtable)
  if (is.null(file)) {
    print(xtable(df, caption = caption, label = label, caption.placement = "top", 
                 align = align, digits = digits),
          sanitize.text.function = function(x){x})
  } else {
    print(xtable(df, caption = caption, label = label, caption.placement = "top", 
                 align = align, digits = digits),
          sanitize.text.function = function(x){x}, file=file, append=TRUE)
  }
}

# hasGroup, specify if values in the last column are groups, it affects how to process matrix
# hasGroup=TRUE, return a data frame by removing last column (groups), 
# and another 1-column data frame for the last column (groups). 
# return a list, where data frame communityMatrix, cols are samples, rows are OTUs/taxa, no preprocessing
# data frame groups may be NULL depending on hasGroup
readCommunityMatrixFile <- function(file, hasGroup) { 
  if(missing(hasGroup)) hasGroup=FALSE

  communityMatrix <- readFile(file)
  cat("\nUpload community matrix file : ", ncol(communityMatrix), "columns,", nrow(communityMatrix), "rows, from", file, "\n") 
  groups <- NULL
  
  # set NA (empty cell) to 0
  communityMatrix[is.na(communityMatrix)] <- 0
  if(hasGroup) {
    groups <- data.frame(row.names=rownames(communityMatrix), groups=communityMatrix[,ncol(communityMatrix)])
    communityMatrix <- communityMatrix[,-ncol(communityMatrix)]
    groups.unique <- unique(groups[,1])
    cat("Split last column to data frame groups that contains", length(groups.unique), 
        "groups : ", paste(groups.unique, collapse=", "), ".\n") 
  } 

  # Return a list 
  list(
    communityMatrix = communityMatrix,
    groups = groups
  )
}

readTaxaFile <- function(file) { 
  taxa <- readFile(file)
  cat("\nUpload taxa file : ", ncol(taxa), "columns,", nrow(taxa), "rows, from", file, "\n") 
  return(taxa)
}

readEnvDataFile <- function(file) { 
  envData <- readFile(file)
  cat("\nUpload environmental data file : ", ncol(envData), "columns,", nrow(envData), "rows, from", file, "\n") 
  return(envData)
}



