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
#' @param in.file The fasta file of OTU representive sequences. 
#' Read by \code{\link{readFasta}}.
#' @param out.file The output fasta file containing extrated sequences.
#' @param otus.names The vector of names to match sequence labels in \code{in.file}.
#' @param regex1,regex2 Use for \code{\link{gsub}(regex1, regex2, id(fasta))} 
#' to remove or replace annotation from original labels. 
#' Default to \code{regex1="(\\|[0-9]+)", regex2=""} in \code{subsetSequences}, 
#' which removes size annotation seperated by "|", 
#' but NA in \code{renameFastaID}, which does nothing.
#' @param ignore.case Refer to \code{\link{gsub}}.
#' @param max.seq Give the number (\code{max.seq}) of seleceted sequences, 
#' if extracted sequences \code{> max.seq}, then choose \code{max.seq} sequences randomly from them. 
#' Defaul to 0 to ignore it.
#' @keywords util
#' @export
#' @examples 
#' subsetSequences(in.file, otus.names, out.file)
#'
#' @rdname utilsSeq
subsetSequences <- function(in.file, out.file, otus.names=c(), regex1="(\\|[0-9]+)", 
                            regex2="", ignore.case = TRUE, max.seq=0) {
  require(ShortRead)
  otus.fasta <- readFasta(in.file)
  if (!is.na(regex1)) {
    # remove/replace annotation
    otus.fasta <- ComMA::renameFastaID(otus.fasta, regex1=regex1, regex2=regex2, ignore.case=ignore.case)
  }
  fasta.id <- as.character(id(otus.fasta))
  
  # extract seq by otus.names
  if (length(otus.names) > 0) {
    if (! is(otus.names, "character"))
      otus.names <- as.character(otus.names)
    
    final.fasta <- otus.fasta[(fasta.id %in% otus.names),]
    
    cat("Extract", length(final.fasta), "sequences from total", length(otus.fasta), 
        "by given", length(otus.names), "sequence names.\n")
    
    if (length(final.fasta) != length(otus.names))
      warning("The extracted sequences ", length(final.fasta), " != ", length(otus.names), " given otus.names !\n")
  } else {
    final.fasta <- otus.fasta
  }
  
  # random select max.seq sequences
  if (max.seq > 1 && length(final.fasta) > max.seq) {
    final.id <- sample(1:length(final.fasta), max.seq, replace=F)
    final.fasta <- final.fasta[final.id]
    cat("Over max limit of sequences ", max.seq, ", so as to randomly extract", length(final.fasta), "sequences from them.\n")
  }
  
  writeFasta(final.fasta, out.file)
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
#' \code{rmDuplicateSeq} removes duplicate sequences or alginments. 
#' 
#' @keywords util
#' @export
#' @examples 
#' rmDuplicateSeq("alg.fasta", "unique-alg.fasta")
#'
#' @rdname utilsSeq
rmDuplicateSeq <- function(in.file, out.file="unique-alg.fasta") {
  require(ShortRead)
  otus.alg <- readFasta(in.file)
  cat("Find", length(srduplicated(otus.alg)[srduplicated(otus.alg)==TRUE]), "duplicate sequences from", in.file, "\n")
  #id(otus.alg)[srduplicated(otus.alg)]
  otus.alg.unique <- otus.alg[srduplicated(otus.alg)==FALSE]
  writeFasta(otus.alg.unique, out.file)
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
getTaxaMap <- function(in.file, trait.name="taxon", trait.value="Bacteria", regex1=NA, regex2="", ignore.case=TRUE) {
  require(ShortRead)
  otus.fasta <- readFasta(in.file)
  # rename the sequences id if regex1 not NA
  otus.fasta <- ComMA::renameFastaID(otus.fasta, regex1=regex1, regex2=regex2, ignore.case=ignore.case)
  
  # 1st col is tip.label
  taxa.map <- data.frame(tip.label=as.character(id(otus.fasta)), stringsAsFactors = F)
  taxa.map[,trait.name] <- trait.value
  return(taxa.map)
}


