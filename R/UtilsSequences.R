# Extract OTUs or reads from CM and taxa classification  
#
# Author: Walter Xie
# Accessed on 20 Apr 2016

#source("http://bioconductor.org/biocLite.R")
#biocLite("ShortRead")
#browseVignettes("ShortRead")

# extract sequences from given otus file whose names match given rownames(cm)   
#' @name utilsSeq
#' @title Utils to manipulate sequences or OTUs
#'
#' @description Utils to manipulate sequences or OTUs, 
#' such as extract OTUs given the subset of OTU names.
#' They are depended on Bioconductor packages  
#' \url{http://www.bioconductor.org/}.
#' 
#' @details 
#' \code{subsetSequences} returns the subset of OTUs matching given OTU names. 
#' It is depended on \code{\link{ShortRead}} package. Follow the instruction 
#' \url{https://bioconductor.org/packages/release/bioc/html/ShortRead.html} to install.
#' 
#' @param otus.file The fasta file of OTU representive sequences. 
#' Read by \code{\link{readFasta}}.
#' @param otus.names The vector of names to match sequence labels in \code{otus.file}.
#' @param otus.out.file The output fasta file containing extrated sequences.
#' @param regex1,regex2 Use for \code{\link{gsub}(regex1, regex2, id(fasta))} 
#' to remove or replace annotation from original labels. 
#' Default to \code{regex1="(\\|[0-9]+)", regex2=""} in \code{subsetSequences}, 
#' which removes size annotation seperated by "|", 
#' but NA in \code{renameFastaID}, which does nothing.
#' @param ignore.case Refer to \code{\link{gsub}}.
#' @param max.seq Give the number (\code{max.seq}) of seleceted sequences randomly in the end. 
#' Defaul to NA to ignore it.
#' @keywords util
#' @export
#' @examples 
#' subsetSequences(otus.file, otus.names, otus.out.file)
#'
#' @rdname utilsSeq
subsetSequences <- function(otus.file, otus.names, otus.out.file, regex1="(\\|[0-9]+)", 
                            regex2="", ignore.case = TRUE, max.seq=NA) {
  require(ShortRead)
  otus.fasta <- readFasta(otus.file)
  # remove/replace annotation
  otus.fasta <- ComMA::renameFastaID(otus.fasta, regex1=regex1, regex2=regex2, ignore.case=ignore.case)
  fasta.id <- as.character(id(otus.fasta))
  
  if (! is(otus.names, "character"))
    otus.names <- as.character(otus.names)
  
  final.fasta <- otus.fasta[(otus.id %in% otus.names),]
  
  cat("Extract", length(final.fasta), "sequences from total", length(otus.fasta), 
      "by given", length(otus.names), "sequence names.\n")
  
  if (length(final.fasta) != length(otus.names))
    stop(cat("The number of extracted sequences", length(final.fasta), 
             "!=", length(otus.names), "given names !\n"))
  
  if (! is.na(max.seq) && max.seq > 1 && length(final.fasta) > max.seq) {
    final.id <- sample(1:length(final.fasta), max.seq, replace=F)
    final.fasta <- final.fasta[final.id]
    cat("Set max sequences = ", max.seq, ", so as to randomly extract", length(final.fasta), "sequences in the end.\n")
  }
  
  writeFasta(final.fasta, otus.out.file)
}

#' @details 
#' \code{renameFastaID} renames the sequences id loaded from a fasta file 
#' using \code{\link{ShortRead}} package. 
#' 
#' @param fasta The fasta object returned by \code{\link{readFasta}}.
#' @keywords util
#' @export
#' @examples 
#' fasta <- renameFastaID(fasta)
#'
#' @rdname utilsSeq
renameFastaID <- function(fasta, regex1=NA, regex2="", ignore.case=TRUE) {
  require(ShortRead)
  if (! is(fasta, "ShortRead") )
    stop("The input has to be ShortRead object !\n", 
         "http://bioconductor.org/packages/release/bioc/html/ShortRead.html")
  
  if (! is.na(regex1)) {
    fasta.id <- gsub(regex1, regex2, id(fasta), ignore.case = ignore.case)
    # change to new id
    fasta.id <- BStringSet(fasta.id)
    fasta@id <- fasta.id
  } 
  return(fasta)
}

#' @details 
#' \code{getTaxaMap} returns a data frame of mapping file to annotate a tree, 
#' such as \code{\link{annotateRAXMLTree}}. 
#' 
#' @param trait.name,trait.value The trait to annotate a tree. 
#' Its format in the string will look like '[trait.name=trait.value]'.
#' @keywords util
#' @export
#' @examples 
#' taxa.map <- getTaxaMap("16s-otus-Bacteria.fasta", trait.name="taxon", trait.value="Bacteria")
#'
#' @rdname utilsSeq
getTaxaMap <- function(otus.file, trait.name="taxon", trait.value="Bacteria", regex1=NA, regex2="", ignore.case=TRUE) {
  require(ShortRead)
  otus.fasta <- readFasta(otus.file)
  # rename the sequences id if regex1 not NA
  otus.fasta <- ComMA::renameFastaID(otus.fasta, regex1=regex1, regex2=regex2, ignore.case=ignore.case)
  
  # 1st col is tip.label
  taxa.map <- data.frame(tip.label=as.character(id(otus.fasta)), stringsAsFactors = F)
  taxa.map[,trait.name] <- trait.value
  return(taxa.map)
}


