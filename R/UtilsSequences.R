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
#' @param regex1,regex2 Use for \code{\link{gsub}(regex1, regex2, id(otus.fasta))} 
#' to remove or replace annotation from original sequence labels. 
#' Default to \code{regex1=NULL, regex2=""}.
#' @param ignore.case Refer to \code{\link{gsub}}.
#' @keywords util
#' @export
#' @examples 
#' extractSequences(otus.file, otus.names, otus.out.file)
#'
#' @rdname utilsSeq
subsetSequences <- function(otus.file, otus.names, otus.out.file, regex1=NULL, regex2="", ignore.case = TRUE) {
  require(ShortRead)
  
  otus.fasta <- readFasta(otus.file)
  otus.id <- id(otus.fasta)
  # remove/replace annotation
  if (! is.null(regex1)) 
    otus.id <- gsub(regex1, regex2, id(otus.fasta), ignore.case = ignore.case)
  
  final.fasta <- otus.fasta[(otus.id %in% otus.names),]
  
  cat("Extract", length(final.fasta), "sequences from total", length(otus.fasta), 
      "by given", length(otus.names), "sequence names.\n")
  
  if (length(final.fasta) != length(otus.names))
    stop(cat("The number of extracted sequences", length(final.fasta), 
             "!=", length(otus.names), "given names !\n"))
  
  writeFasta(final.fasta, otus.out.file)
}

