# Author: Walter Xie
# Accessed on 10 Mar 2016

#' Create folder from \code{\link{file.path}} if not exist.
#' 
#' @param dir.path The folder path constructed by \code{\link{file.path}}.
#' @keywords IO
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
split_path <- function(path) {
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
#' communityMatrix <- readFile("16S.txt", msg.file="16S OTU table", msg.col="samples", msg.row="OTUs")
#' taxa.table <- readFile("16S_taxonomy_table.txt", msg.file="16S taxonomy table", msg.row="OTUs")
#' env <- readFile("env_data.txt", msg.file="enviornmental data", msg.row="samples")
#' taxa.phyla <- readFile("taxonomy97phyla.txt", row.names=NULL)
readFile <- function(file, sep="\t", header=TRUE, row.names=1, verbose=TRUE, msg.file="file", msg.col="columns", msg.row="rows") { 
  require(tools)
  # sep="\t" only work for non csv file
  if (tolower(file_ext(file))=="csv") {
    df <- read.csv(file, header=header, row.names=row.names, check.names=FALSE, stringsAsFactors=FALSE)
  } else {
    # be careful read.table bug   
    df <- read.table(file, sep=sep, header=header, row.names=row.names, check.names=FALSE, stringsAsFactors=FALSE)  
  }	
  
  if (verbose) 
    cat("\nUpload", msg.file, ":", ncol(df), msg.col, ",", nrow(df), msg.row, ", from", file, "\n") 
  
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
                       na = "NA", dec = ".", row.names = TRUE, col.names = TRUE) {
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
#' @keywords IO
#' @examples 
#' tableFile <- file.path(workingPath, "report.tex")
#' printXTable(data.frame, caption = "Phylogenetic beta diversity", 
#'             label = "tab:pd:beta", file=tableFile)
printXTable <- function(df, caption, label, include.rownames=TRUE, file=NULL, align = NULL, digits = NULL) {
  require(xtable)
  if (is.null(file)) {
    print(xtable(df, caption = caption, label = label, caption.placement = "top", 
                 align = align, digits = digits),
          sanitize.text.function = function(x){x}, include.rownames=include.rownames)
  } else {
    print(xtable(df, caption = caption, label = label, caption.placement = "top", 
                 align = align, digits = digits),
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

# hasGroup, specify if values in the last column are groups, it affects how to process matrix
# hasGroup=TRUE, return a data frame by removing last column (groups), 
# and another 1-column data frame for the last column (groups). 
# return a list, where data frame communityMatrix, cols are samples, rows are OTUs/taxa, no preprocessing
# data frame groups may be NULL depending on hasGroup
#' @name ComMAIO
#' @title IO functions in \pkg{ComMA}
#'
#' @description IO functions for the data defined in \pkg{ComMA} package, such as community matrix.
#' 
#' @details 
#' \code{readCommunityMatrix} reads file to return a community matrix.
#' 
#' @param file The file to read.
#' @keywords IO
#' @export
#' @examples 
#' cm <- readCommunityMatrix("16S.txt", "16S")
#' 
#' @rdname ComMAIO
readCommunityMatrix <- function(file, matrix.name=NULL, minAbund=2, regex="(\\|[0-9]+)", verbose=TRUE) { 
  communityMatrix <- ComMA::readFile(file, verbose=verbose, msg.file=paste(matrix.name, "community matrix"), 
                                     msg.col="samples", msg.row="OTUs")
  communityMatrix <- ComMA::rmMinAbundance(communityMatrix, minAbund)
  
  if (!is.null(regex))
    rownames(communityMatrix) <- gsub(regex, "", rownames(communityMatrix))
  
  attr(communityMatrix,"name") <- matrix.name
  return(communityMatrix)
}

# "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"
# "kingdom" is compulsory, "path" is optional
# c("16s-path.txt", "16s-kingdom.txt", "16s-phylum.txt", "16s-class.txt", "16s-order.txt", "16s-family.txt")
#' \code{taxaTableMEGAN} reads a list of taxa mapping files exported from MEGAN to create a taxa table in \code{folder.path}.
#' 
#' @param file.prefix The prefix to find taxa mapping files, 
#' such as "16s" for "16s-path.txt", "16s-kingdom.txt", "16s-phylum.txt".
#' @param folder.path The folder path to contain taxa mapping files 
#' and output MEGAN taxa table.
#' @keywords IO
#' @export
#' @examples 
#' taxaTableMEGAN("16s", getwd())
#' 
#' @rdname ComMAIO
taxaTableMEGAN <- function(file.prefix, folder.path, col.names=c("path", "kingdom", "phylum", "class", "order", "family"), 
                           sep="\t", regex="(\\|[0-9]+)") {
  taxa.files <- paste0("16s-", col.names, ".txt")
  
  for (n in 1:length(col.names)) {
    t.f <- file.path(folder.path, taxa.files[n])
    # no header, no row.names in the file
    df.taxa <- readFile(t.f, sep=sep, header=FALSE, row.names=NULL) 
    if (ncol(df.taxa) < 2)
      stop(paste("Taxa file", t.f, "can be correctly parsed ! Please check the file format."))
    
    colnames(df.taxa) <- c("OTUs",col.names[n])
    if (n==1)
      df.path <- df.taxa
    else 
      df.path <- merge(df.path, df.taxa, by = "OTUs")
  }
  
  if (!is.null(regex))
    df.path[,"OTUs"] <- gsub(regex, "", df.path[,"OTUs"])
  
  megan.f <- file.path(folder.path, "16s-megan.txt")
  cat("Write taxa table", nrow(df.path), "rows", ncol(df.path), "columns to file", megan.f,".\n")
  writeTable(df.path, megan.f, row.names = FALSE)
}

# Zxan08_H415I8K02GEDIZ|2	k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__;g__;s__	0.970
taxaTableRDP <- function(file, minCred=0.8, sep="\t", regex="(\\|[0-9]+)") {
  df.taxa <- readFile(file, sep=sep, header=FALSE, row.names=NULL) 
  
  if (ncol(df.taxa) < 3)
    stop(paste("RDP output file", file, "can be correctly parsed ! Please check the file format."))
  
  colnames(df.taxa) <- c("OTUs","path","credibility")
  # assign unclassified
  df.taxa[df.taxa[,"credibility"] < minCred, "path"] <- "unclassified"
  cat("Assign", nrow(df.taxa[df.taxa[,"path"] == "unclassified",]), 
      "rows to 'unclassified' whose credibility <", minCred, ".\n")
  
  
  for (n in 1:length(taxaFiles)) {
    
    colnames(df.taxa) <- c("OTUs",ranks[n])
    if (n==1)
      df.path <- df.taxa
    else 
      df.path <- merge(df.path, df.taxa, by = "OTUs")
  }
  if (!is.null(regex))
    df.path[,"OTUs"] <- gsub(regex, "", df.path[,"OTUs"])
  
  cat("Write taxa table", nrow(df.path), "rows", ncol(df.path), "columns to file", ,".\n")
  
  writeTable(df.path, file.path(folder.path, "16s-megan.txt"), row.names = FALSE)
}


#' \code{readTaxaTable} reads a file to return a taxa table
#' 
#' @return 
#' A taxa table.
#' @keywords IO
#' @export
#' @examples 
#' tt.megan <- readTaxaTable("16s-megan.txt", "16S taxa table", taxa.group="Bacteria")
#' 
#' @rdname ComMAIO
readTaxaTable <- function(file, matrix.name=NULL, taxa.group="assigned", rank="kingdom", include=TRUE, verbose=TRUE) {
  taxa.table <- ComMA::readFile(file, verbose=verbose,  msg.file=paste(matrix.name, "taxonomy table"), msg.row="OTUs")
  taxa.table <- taxa.table[order(rownames(taxa.table)),]
  # make lower case to match ranks
  colnames(taxa.table) <- tolower(colnames(taxa.table))
  
  n.taxa=nrow(taxa.table)
  ##### keep OTU rows contain given taxa belongTo ##### 
  if (taxa.group != "all") {
    # Exclude unassigned etc
    taxa.table <- subset(taxa.table, !(grepl("root|cellular organisms|No hits|Not assigned", taxa.table[,"kingdom"], ignore.case = T)))  

    if (taxa.group != "assigned") {
      if (include) {
        # for PROTISTS, taxa.group="CHROMISTA|PROTOZOA", rank="kingdom"
        taxa.table <- subset(taxa.table, (grepl(taxa.group, taxa.table[,rank], ignore.case = T))) 
      } else { 
        # exclude some phyla, taxa.group="Cnidaria|Brachiopoda|Echinodermata|Porifera", rank="phylum"
        taxa.table <- subset(taxa.table, !grepl(taxa.group, taxa.table[,rank], ignore.case = T)) 
      }
    }
  }
  
  if (nrow(taxa.table) < 1)
    cat("Warning: cannot find", taxa.group, "at", rank, "from taxa path file", file, "!")
  
  if(verbose && nrow(taxa.table) < n.taxa) 
    cat("\nSelect", nrow(taxa.table), "classifications, taxa.group =", taxa.group, 
        ", rank =", rank, ", include =", include, ".\n") 
  
  return(taxa.table)
}

