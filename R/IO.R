# Author: Walter Xie
# Accessed on 10 Mar 2016

#' Create folder from \code{\link{file.path}} if not exist.
#' 
#' @param dir.path The folder path constructed by \code{\link{file.path}}.
#' @keywords IO
#' @export
#' @examples 
#' figDir <- file.path(workingPath, "figures")
#' ComMA::mkdir(figDir) 
mkDir <- function(dir.path) {
  require(tools)
  if (!file.exists(dir.path)) {    
    dir.create(dir.path)    
  }
  cat("\nConfig : setup", dir.path, "\n")
}

#' Split a file path into folder names vector
#' 
#' @param path The file path to be split.
#' @source \url{http://stackoverflow.com/questions/29214932/split-a-file-path-into-folder-names-vector}
#' @return 
#' A the reversed vector of folder names.
#' @keywords IO
#' @export
#' @examples 
#' split_path("/home/foo/stats/index.html")
#' [1] "index.html" "stats"      "foo"        "home"      
#' split_path("C:\\Windows\\System32")
#' [1] "System32" "Windows"  "C:"      
splitPath <- function(path) {
  rev(setdiff(strsplit(path,"/|\\\\")[[1]], ""))
}

#' Read a file to a data frame. 
#' 
#' If the file extension is \emph{csv}, 
#' then use \code{\link{read.csv}}, otherwise use \code{\link{read.table}}.
#' 
#' @param file The 1st row is column names, the 1st column is row names.
#' @param row.names A vector of row names. This can be a vector giving the actual row names, 
#' or a single number giving the column of the table which contains the row names, 
#' or character string giving the name of the table column containing the row names.
#' If there is a header and the first row contains one fewer field than the number of columns, 
#' the first column in the input is used for the row names. 
#' Otherwise if row.names is NULL, the rows are numbered. Default to 1.
#' Using row.names = NULL forces row numbering.
#' @param header A logical value indicating whether the file contains 
#' the names of the variables as its first line. Default to TRUE. 
#' @param sep Only used for non \emph{csv} file. Default to tab "\\t".
#' @param msg.file,msg.col,msg.row The message regarding file, column and row, if verbose=TRUE.  
#' @param verbose More details. Default to TRUE.
#' @return 
#' A data frame from the file, such as
#' \tabular{rrrr}{
#'   OTU_id \tab plot01 \tab plot02\tab ...\cr
#'   OTU_1 \tab 1 \tab 0 \tab ...\cr
#'   OTU_2 \tab 100 \tab 200 \tab ...\cr
#'   OTU_3 \tab 56 \tab 3 \tab ...
#' }
#' @export
#' @keywords IO
#' @examples 
#' community.matrix <- readFile("16S.txt", msg.file="16S OTU table", msg.col="samples", msg.row="OTUs")
#' taxa.table <- readFile("16S_taxonomy_table.txt", msg.file="16S taxonomy table", msg.row="OTUs")
#' env <- readFile("env_data.txt", msg.file="enviornmental data", msg.row="samples")
#' taxa.phyla <- readFile("taxonomy97phyla.txt", row.names=NULL)
readFile <- function(file, sep="\t", header=TRUE, row.names=1, verbose=TRUE, 
                     msg.file="file", msg.col="columns", msg.row="rows") { 
  require(tools)
  # sep="\t" only work for non csv file
  if (tolower(file_ext(file))=="csv") {
    df <- read.csv(file, header=header, row.names=row.names, check.names=FALSE, stringsAsFactors=FALSE)
  } else {
    # be careful read.table bug   
    df <- read.table(file, sep=sep, header=header, row.names=row.names, check.names=FALSE, stringsAsFactors=FALSE)  
  }	
  
  if (verbose) {
    if (header) 
      msg.row <- paste(msg.row, "excluding column names from row 1")
    if (!is.null(row.names))
      msg.col <- paste(msg.col, "excluding row names from column", row.names)
    
    cat("\nUpload", msg.file, ":", ncol(df), paste0(msg.col, ","), nrow(df), paste0(msg.row, ","), "from file", file, "\n") 
  }
  
  return(df)
}

#' Write the data fram (table) to a file.
#' 
#' @param df A data frame. 
#' @param file If the file extension is \emph{csv}, 
#' then use \code{\link{write.csv}}, 
#' otherwise use \code{\link{write.table}}.
#' @param ... More parameters, see \code{\link{write.table}}.
#' @export
#' @keywords IO
#' @examples 
#' writeTable(df, file.path)
writeTable <- function(df, file, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", 
                       na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, 
                       verbose=TRUE, msg.file="file", msg.col="columns", msg.row="rows") {
  require(tools)
  if (tolower(file_ext(file))=="csv") {
    write.csv(df, file, append = append, quote = quote, eol = eol, 
              na = na, dec = dec, row.names = row.names, col.names = col.names)
  } else if (row.names == TRUE && col.names == TRUE) {
    cat(paste(c("row.names", colnames(df)), collapse = sep), file = file, sep = "\n")
    write.table(df, file, append = TRUE, quote = quote, eol = eol, sep = sep,
                na = na, dec = dec, row.names = TRUE, col.names = FALSE)
  } else {# .tsv .txt
    #write.table bug: if row.names = TRUE, mistake to start col names from the 1st cell
    write.table(df, file, append = append, quote = quote, eol = eol, sep = sep,
                na = na, dec = dec, row.names = row.names, col.names = col.names)
  }
  
  if (verbose) {
    cat("\nWrite", msg.file, ":", ncol(df), paste0(msg.col, ","), nrow(df), paste0(msg.row, ","), "to file", file, "\n") 
  }
}

#' Print data fram (table) as Latex format to either file or console, 
#' which is extened from \code{\link{xtable}}.
#' 
#' @param df A data frame.
#' @param caption Latex table caption.
#' @param label Latex table label.
#' @param file If NULL, then print the results to console, 
#' otherwise print them to the file. Default to NULL. 
#' @param invalid.char If it is TRUE, then add \\\\ to all 
#' _ or % for Latex. Default to TRUE.
#' @param include.rownames,caption.placement 
#' Refer to \code{\link{print.xtable}}.
#' @param ... More at \code{\link{xtable}}. 
#' @export
#' @keywords IO
#' @examples 
#' tableFile <- file.path(workingPath, "report.tex")
#' printXTable(data.frame, caption = "Phylogenetic beta diversity", 
#'             label = "tab:pd:beta", file=tableFile)
printXTable <- function(df, caption="The table of", label="tab:label", file=NULL, 
                        invalid.char=TRUE, include.rownames=TRUE, caption.placement="top", ...) {
  if (invalid.char) {
    rownames(df) <- gsub("_", "\\\\_", rownames(df))
    colnames(df) <- gsub("_", "\\\\_", colnames(df))
    df <- as.data.frame(gsub("%|_", "\\\\%|\\\\_", as.matrix(df)))
  }
    
  require(xtable)
  if (is.null(file)) {
    print(xtable(df, caption = caption, label = label, caption.placement = caption.placement, ...),
          sanitize.text.function = function(x){x}, include.rownames=include.rownames)
  } else {
    print(xtable(df, caption = caption, label = label, caption.placement = caption.placement, ...),
          sanitize.text.function = function(x){x}, include.rownames=include.rownames, 
          file=file, append=TRUE)
  }
}

# for Latex file
mkValidTex <- function(file) {
  # change all _, %
  system(paste("sed -i.bak 's/\\_/\\\\_/g;' ", report.file))
  system(paste("sed -i.bak 's/\\%/\\\\%/g;' ", file)) 
}


#' Read File Line By Line
#' 
#' @param file The file to read.
#' @return 
#' A character vector of length the number of lines read. Refer to \code{\link{readLines}}.
#' @keywords IO
#' @export
#' @examples 
#' linn <- readFileLineByLine("time.txt")
readFileLineByLine <- function(file) {
  conn <- file(file,open="r")
  linn <-readLines(conn)
  close(conn)
  return(linn)
}
