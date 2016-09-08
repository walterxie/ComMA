# Author: Walter Xie
# Accessed on 10 Mar 2016


# hasGroup, specify if values in the last column are groups, it affects how to process matrix
# hasGroup=TRUE, return a data frame by removing last column (groups), 
# and another 1-column data frame for the last column (groups). 
# return a list, where data frame community.matrix, cols are samples, rows are OTUs/taxa, no preprocessing
# data frame groups may be NULL depending on hasGroup
#' @name IOComMA
#' @title IO functions in \pkg{ComMA}
#'
#' @description IO functions for the data defined in \pkg{ComMA} package, such as community matrix.
#' 
#' @details 
#' \code{readCommunityMatrix} reads file to return a community matrix.
#' 
#' @param file The file to read/write.
#' @param matrix.name The string to locate the matrix from its file name.
#' @param minAbund The minimum abundance threshold to remove rows/columns 
#' by row/column sum of abundance. For exampe, if minAbund=2, then remove 
#' all singletons appeared in only one sample. If minAbund=1, 
#' then remove all empty rows/columns. Default to 2 (singletons).
#' But \code{postfix} is only used for naming, no data processed.
#' @param regex1,regex2 Use for \code{\link{gsub}(regex1, regex2, row.names)} 
#' to remove or replace annotation from original labels. 
#' Default to \code{regex1="(\\|[0-9]+)", regex2=""} for read??? but NULL to write???, 
#' which removes size annotation seperated by "|".
#' @param ignore.case Refer to \code{\link{gsub}}.
#' @param col.name.decreasing Should the sort decreasing order of \code{colnames} 
#' be TRUE? Refer to \code{\link{order}}. If NULL, do nothing.
#' Default to FALSE.
#' @keywords IO
#' @export
#' @examples 
#' cm <- readCommunityMatrix("16S.txt", "16S")
#' 
#' @rdname IOComMA
readCommunityMatrix <- function(file, matrix.name=NULL, minAbund=2, 
                                regex1="(\\|[0-9]+)", regex2="", ignore.case=TRUE,
                                col.name.decreasing=FALSE, verbose=TRUE) { 
  community.matrix <- ComMA::readFile(file, verbose=verbose, msg.file=paste(matrix.name, "community matrix"), 
                                     msg.col="samples", msg.row="OTUs")
  community.matrix <- ComMA::rmMinAbundance(community.matrix, minAbund)
  
  # remove/replace annotation
  if (! is.null(regex1)) 
    rownames(community.matrix) <- gsub(regex1, regex2, rownames(community.matrix), ignore.case = ignore.case)
  
  if (! is.null(col.name.decreasing)) {
    community.matrix <- community.matrix[,order(colnames(community.matrix), decreasing=col.name.decreasing)]   
  } 
  
  attr(community.matrix,"name") <- matrix.name
  return(community.matrix)
}

#' @details 
#' \code{writeCommunityMatrix} writes a community matrix to file.
#' 
#' @keywords IO
#' @export
#' @examples 
#' writeCommunityMatrix(cm, "16S-new.txt", msg.file="16S")
#' 
#' @rdname IOComMA
writeCommunityMatrix <- function(community.matrix, file, 
                                 msg.file="file", msg.col="columns", msg.row="rows",
                                 regex1=NULL, regex2="", ignore.case=TRUE) { 
  # rm empty rows and columns
  community.matrix <- ComMA::rmMinAbundance(community.matrix, minAbund=1)
  community.matrix <- ComMA::rmMinAbundance(community.matrix, minAbund=1, MARGIN=2)
  
  # remove/replace annotation
  if (! is.null(regex1)) 
    rownames(community.matrix) <- gsub(regex1, regex2, rownames(community.matrix), 
                                       ignore.case = ignore.case)
  
  ComMA::writeTable(community.matrix, file, msg.file=msg.file, msg.col=msg.col, msg.row=msg.row)
}

#' @details 
#' \code{readTaxaTable} reads a file to return a taxa table. 
#' \code{\link{subsetTaxaTable}} can be used to takes or excludes 
#' a subset of the taxa table recursively.
#' 
#' @param taxa.group The taxonomic group, the values can be 'all', 'assigned', or 
#' Group 'all' includes everything.
#' Group 'assigned' removes all uncertain classifications including 
#' 'root', 'cellular organisms', 'No hits', 'Not assigned'. 
#' Alternatively, any high-ranking taxonomy in your taxonomy file 
#' can be used as a group or multi-groups (seperated by "|"), 
#' such as 'BACTERIA', 'Proteobacteria', etc. But they have to be 
#' in the same rank column in the file. Default to remove all 
#' uncertain classifications, even when group(s) assigned.
#' @param rank The rank column in the file to specify where to 
#' search for \code{taxa.group}. 
#' @param include Define whether include or exclude given \code{taxa.group}. 
#' Default to TRUE.
#' @keywords IO
#' @export
#' @examples 
#' tt.megan <- readTaxaTable("16s-megan.txt", "16S taxa table", taxa.group="Bacteria")
#' tt.rdp.8 <- readTaxaTable("16s-rdp.txt", "16S taxa table", taxa.group="Bacteria")
#' 
#' @rdname IOComMA
readTaxaTable <- function(file, matrix.name=NULL, taxa.group="assigned", rank="kingdom", 
                          regex1="(\\|[0-9]+)", regex2="", ignore.case=TRUE, 
                          include=TRUE, sep="\t", verbose=TRUE) {
  taxa.table <- ComMA::readFile(file, sep=sep, verbose=verbose, 
                                msg.file=paste(matrix.name, "taxonomy table"), msg.row="OTUs")
  taxa.table <- taxa.table[order(rownames(taxa.table)),]
  # make lower case to match ranks
  colnames(taxa.table) <- tolower(colnames(taxa.table))
  
  if (! is.element(rank, colnames(taxa.table)) )
    stop("Cannot find rank column ", rank, "! Use createTaxaTable??? function to create taxa.table file.")
  
  # remove/replace annotation
  if (! is.null(regex1)) 
    rownames(taxa.table) <- gsub(regex1, regex2, rownames(taxa.table), ignore.case = ignore.case)
  
  n.taxa=nrow(taxa.table)
  ##### keep OTU rows contain given taxa belongTo ##### 
  if (taxa.group != "all") {
    # Exclude unassigned etc
    taxa.table <- subset(taxa.table, !(grepl("root|cellular organisms|No hits|Not assigned", taxa.table[,"kingdom"], ignore.case = T)))  
    
    # for PROTISTS, taxa.group="CHROMISTA|PROTOZOA"
    if (toupper(taxa.group) == "PROTISTS")
      taxa.group="CHROMISTA|PROTOZOA"
    
    if (taxa.group != "assigned") {
      taxa.table <- ComMA::subsetTaxaTabl(taxa.table, taxa.group, rank, include)
    }
  }
  
  if (nrow(taxa.table) < 1)
    warning("Cannot find ", taxa.group, " at ", rank, " from taxa path file ", file, " !")
  
  if(verbose && nrow(taxa.table) < n.taxa) {
    cat("Select", nrow(taxa.table), "classifications, by given taxa.group =", taxa.group, 
        ", rank =", rank, ", include =", include, ".\n") 
  }
  
  return(taxa.table)
}

# "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"
# "kingdom" is compulsory, "path" is optional
# c("16s-path.txt", "16s-kingdom.txt", "16s-phylum.txt", "16s-class.txt", "16s-order.txt", "16s-family.txt")
#' @details 
#' \code{createTaxaTableMEGAN} reads a list of taxa mapping files exported 
#' from MEGAN to create a taxa table in \code{folder.path}.
#' 
#' @param file.prefix The prefix to find taxa mapping files, 
#' such as "16s" for "16s-path.txt", "16s-kingdom.txt", "16s-phylum.txt".
#' @param folder.path The folder path to contain taxa mapping files 
#' and output MEGAN taxa table.
#' @param col.names A vector of column names of taxonomic ranks in the taxa table, 
#' and they also determine the files to be merged into columns. Default to 
#' c("path", "kingdom", "phylum", "class", "order", "family", "genus"), which indicates 
#' the list of files c("16s-path.txt", "16s-kingdom.txt", "16s-phylum.txt", "16s-class.txt", 
#' "16s-order.txt", "16s-family.txt", "16s-genus.txt") given file.prefix = "16s". "path" is optional.
#' @param file.out The taxonomic table file name.
#' @keywords IO
#' @export
#' @examples 
#' createTaxaTableMEGAN(getwd(), "16s", file.out="16s-megan.txt")
#' 
#' @rdname IOComMA
createTaxaTableMEGAN <- function(folder.path, file.prefix, sep="\t", 
                                 regex1="(\\|[0-9]+)", regex2="", ignore.case=TRUE, 
                                 col.names=c("path", "kingdom", "phylum", "class", "order", "family", "genus"), 
                                 file.out="taxa-table-megan.txt") {
  taxa.files <- paste0(file.prefix, "-", col.names, ".txt")
  
  for (n in 1:length(col.names)) {
    t.f <- file.path(folder.path, taxa.files[n])
    # no header, no row.names in the file
    df.taxa <- ComMA::readFile(t.f, sep=sep, header=FALSE, row.names=NULL) 
    if (ncol(df.taxa) < 2)
      stop(paste("Taxa file", t.f, "can be correctly parsed ! Please check the file format."))
    
    # remove " <phylumn>" added from MEGAN
    df.taxa[,2] <- gsub("\\s+<.*>", "", df.taxa[,2])
    # remove " (miscellaneous)"
    df.taxa[,2] <- gsub("\\s+\\(miscellaneous\\)", "", df.taxa[,2])

    colnames(df.taxa) <- c("OTUs",col.names[n])
    if (n==1)
      df.path <- df.taxa
    else 
      df.path <- merge(df.path, df.taxa, by = "OTUs")
  }
  
  # remove/replace annotation
  if (! is.null(regex1)) 
    df.path[,"OTUs"] <- gsub(regex1, regex2, df.path[,"OTUs"], ignore.case = ignore.case)
  
  taxa.f <- file.path(folder.path, file.out)
  ComMA::writeTable(df.path, taxa.f, row.names = FALSE, msg.file = "taxa table")
}

#' @details 
#' \code{createTaxaTableRDP} reads a 3-column taxa mapping file generated from RDP 
#' to create a taxa table in \code{folder.path}. The row look like: 
#' OTU1	k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__;g__;s__	0.970
#' 
#' @param rm.rank.prefix Remove rank prefix, such as 'k__'. Default to TRUE.
#' @keywords IO
#' @export
#' @examples 
#' createTaxaTableRDP("otus_tax_assignments.txt", file.out="16s-rdp.txt")
#' 
#' @rdname IOComMA
createTaxaTableRDP <- function(file, sep="\t", rm.rank.prefix=TRUE, 
                               regex1="(\\|[0-9]+)", regex2="", ignore.case=TRUE,  
                               file.out="taxa-table-rdp.txt") {
  df.taxa <- ComMA::readFile(file, sep=sep, header=FALSE, row.names=NULL) 
  if (ncol(df.taxa) < 3)
    stop(cat("RDP output file", file, "can be correctly parsed !\nPlease check the file format."))
  
  colnames(df.taxa) <- c("OTUs","path","confidence")
  # rm []
  df.taxa[,"path"] <- gsub("\\[|\\]", "", df.taxa[,"path"])
  # rm Root; for early version 
  df.taxa[,"path"] <- gsub("Root;", "", df.taxa[,"path"], ignore.case = T)
  
  if (rm.rank.prefix) 
    vector.path <- gsub("[a-z]__", "", df.taxa[,"path"])
  else 
    vector.path <- df.taxa[,"path"] 
  
  cols.lin <- read.table(text = vector.path, sep = ";", colClasses = "character", fill=TRUE)
  if (ncol(cols.lin) < 3)
    stop(cat("RDP output file", file, ", 2nd column 'taxa path' needs at least 3 ranks !",
             "\nFor example, k__;p__;c__;o__;f__;g__;s__"))
  
  # fix to k__;p__;c__;o__;f__;g__;s__
  colnames(cols.lin) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")[1:ncol(cols.lin)]
  
  df.path <- cbind(df.taxa, cols.lin)

  # remove/replace annotation
  if (! is.null(regex1)) 
    df.path[,"OTUs"] <- gsub(regex1, regex2, df.path[,"OTUs"], ignore.case = ignore.case)
  
  ComMA::writeTable(df.path, file.out, row.names = FALSE, msg.file = "taxa table")
}

