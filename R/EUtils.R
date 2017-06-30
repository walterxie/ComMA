# Author: Walter Xie
# Accessed on 22 Jun 2017

#' @name EUtils
#' @title EUtils NCBI
#'
#' @description Utils to process data downloaded from NCBI.
#' @import reutils
#' @import xml2
#' 
#' @details 
#' \code{parseTaxaSet} parses the taxonomy XML (DOCTYPE is TaxaSet) 
#' as the result of \code{\link{efetch}} from taxonomy database into a data.frame, 
#' which inlcudes "TaxId", "ScientificName", "Rank", "Lineage", "Division", 
#' and the format of taxa.table from "kingdom" to "genus".
#' 
#' <TaxaSet>
#' <Taxon>
#'   <TaxId>123685</TaxId>
#'   <ScientificName>Oryzias minutillus</ScientificName>
#'   <ParentTaxId>8089</ParentTaxId>
#'   <Rank>species</Rank>
#'   <Division>Vertebrates</Division>
#'   <Lineage>cellular organisms; Eukaryota; Opisthokonta; Metazoa; ...</Lineage>
#'   <LineageEx>
#'    <Taxon>
#'     <TaxId>131567</TaxId>
#'     <ScientificName>cellular organisms</ScientificName>
#'     <Rank>no rank</Rank>
#'    </Taxon>
#'    ...  
#'   </LineageEx>
#' </Taxon>
#' ...
#' </TaxaSet>
#' 
#' @param xml The XML result of \code{\link{efetch}}. 
#' @keywords eutils
#' @export
#' @examples 
#' library("reutils")
#' taxa <- efetch(c("123685", "8089", "8088"), "taxonomy")
#' taxa.df <- parseTaxaSet(taxa$content)
#' 
#' @rdname EUtils
parseTaxaSet <- function(xml) {
  tag1 <- "TaxaSet"
  tags <- c("TaxId", "ScientificName", "Rank", "Lineage", "Division")
  ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")
  
  require(xml2)
  x <- read_xml(xml)
  # parse TaxaSet
  taxonset <- xml_children(x)
  if (xml_name(x) != tag1 || length(taxonset) < 1) 
    stop("No child below tag <", tag1, "> !")
  
  taxa.df <- NULL
  # parse each Taxon
  for (t in 1:length(taxonset)) {
    taxon <- xml_children(taxonset[t])
    # parse lineage
    lin.x <- xml_children(taxon)
    lin.taxid <- xml_text(xml_find_all(lin.x, paste0(".//", tags[1]))) # TaxId
    lin.name <- xml_text(xml_find_all(lin.x, paste0(".//", tags[2]))) # ScientificName
    lin.rank <- xml_text(xml_find_all(lin.x, paste0(".//", tags[3]))) # Rank
    # convert to taxa table format
    n.r <- lin.name[match(ranks, lin.rank)]
    
    # insert a row with taxa table format appended in the end
    rbind(taxa.df, c(xml_text(taxon), n.r)) -> taxa.df
    if (t==1)
      colnames(taxa.df) <- c(xml_name(taxon), ranks)
  }
  
  if (!all(tags %in% colnames(taxa.df)))
    stop("Cannot find tags : ", paste(tags, collapse = ","), " !")
  taxa.df <- taxa.df[,c(tags, ranks)]
  data.frame(taxa.df, stringsAsFactors=FALSE)
}

#' @details 
#' \code{parseTSeqSet} parses the TinySeq XML (DOCTYPE is TSeqSet) 
#' as the result of \code{\link{efetch}} from nuccore database into a data.frame, 
#' which inlcudes "TaxId", "ScientificName", "ACCESSION", "Lineage", "sequence".
#' 
#' <TSeqSet>
#'   <TSeq>
#'     <TSeq_seqtype value="nucleotide"/>
#'     <TSeq_gi>1079489517</TSeq_gi>
#'     <TSeq_accver>NC_031445.1</TSeq_accver>
#'     <TSeq_sid>gnl|NCBI_GENOMES|60824</TSeq_sid>
#'     <TSeq_taxid>126358</TSeq_taxid>
#'     <TSeq_orgname>Abeliophyllum distichum</TSeq_orgname>
#'     <TSeq_defline>Abeliophyllum distichum chloroplast, complete genome</TSeq_defline>
#'     <TSeq_length>155982</TSeq_length>
#'     <TSeq_sequence>CATTTTAGTTATGGGC...GCTGT</TSeq_sequence>
#'   </TSeq>
#' </TSeqSet>
#' 
#' @param save.seq If TRUE as default, save sequences 
#' into data.frame returned by \code{parseTSeqSet}.
#' @export
#' @examples 
#' seqset <- efetch("NC_031445.1", "nuccore", "fasta")
#' seqset.df <- parseTSeqSet(seqset$content)
#' 
#' @rdname EUtils
parseTSeqSet <- function(xml, save.seq=T) {
  tag1 <- "TSeqSet"
  
  require(xml2)
  x <- read_xml(xml)
  # parse TSeqSet
  seqset <- xml_children(x)
  if (xml_name(x) != tag1 || length(seqset) < 1) 
    stop("No child below tag <", tag1, "> !")
  
  seq.df <- NULL
  # parse each TSeq
  for (t in 1:length(seqset)) {
    seq <- xml_children(seqset[t])
    
    # insert a row 
    rbind(seq.df, xml_text(seq)) -> seq.df
    if (t==1)
      colnames(seq.df) <- gsub("TSeq_", "", xml_name(seq)) # rm TSeq_ prefix
    if ("seqtype" %in% colnames(seq.df)) {
      seqtype <- xml_child(seqset[t], "TSeq_seqtype")
      seq.df[t,"seqtype"] <- xml_attr(seqtype, "value")
    }
    if (!save.seq && "sequence" %in% colnames(seq.df))
      seq.df[t,"sequence"] <- NA
  }
  
  data.frame(seq.df, stringsAsFactors=FALSE)
}

#' @details 
#' \code{batchDownload} downloads NCBI data given a vector of \code{uid} 
#' using \code{\link{efetch}} one at a time for each \code{uid}, 
#' and writes all results directly to a file without parsing.  
#' 
#' @param uid.vec The vector of \code{uid} for \code{\link{efetch}} to NCBI data. 
#' The \code{uid} will either look like NC_016668.1 ("accver") or KM462867 ("INSDC").
#' @param db,rettype,retmode,... The arguments of \code{\link{efetch}}.
#' See \url{https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly}
#' for all supported databases and their available values. 
#' @param out.file The file to write all results directly without parsing.
#' Note: if \code{length(uid) > 500}, \code{outfile} of \code{\link{efetch}} is required.
#' @param sleep Please be nice to give enough time to break. Default 10 seconds.
#' @keywords eutils
#' @export
#' @examples
#' batchDownload(c("NC_031445.1", "NC_026892.1"), "nuccore", "gb", out.file="res.gbff") 
#'  
#' @rdname EUtils
batchDownload <- function(uid.vec, db = NULL, rettype = NULL, retmode = NULL, out.file="res.txt", sleep=10, ...) {
  if (length(uid.vec) < 0)
    stop("Invalid input vector of uid ! length = ", length(uid.vec))
  
  require("reutils")
  res.vec <- c()
  for (i in 1:length(uid.vec)) {
    uid <- as.character(uid.vec[i])
    cat(i, ": downloading ", uid, " from db =", db, ", retrieval type =", rettype, " \n")
    res <- efetch(uid, db=db, rettype=rettype, retmode=retmode, ...)
    res.vec <- c(res.vec, res$content)
    Sys.sleep(10)
  }
  cat("Write to file ", out.file, "\n")
  file.create(out.file)
  fileConn <- file(out.file)
  writeLines(gb.vec, fileConn)
  close(fileConn)
}

#' @details 
#' \code{seqSet2Fasta} downloads reference sequences given their \code{uid} 
#' using \code{\link{efetch}}. \code{parseTSeqSet} is used to parse the 
#' result of \code{\link{efetch}}. 
#' 
#' @param seqset.uid The vector of \code{uid} for \code{\link{efetch}} to download sequence,
#' where db = "nuccore", rettype = "fasta". 
#' The \code{uid} will either look like NC_016668.1 ("accver") or KM462867 ("INSDC").
#' @param by Split \code{seqset.uid} into subgroups by the given number, defaul to 20.
#' @param folder,file.prefix,file.extension Determine the file name.
#' @param seq.label The vector of string to determine how to label the sequence. 
#' If the element is one of the column name data.frame from \code{parseTSeqSet}, 
#' then the label of that position will be the value of that column. 
#' The available columns are "seqtype", "gi", "accver", "taxid", "orgname", "defline", "length". 
#' @keywords eutils
#' @export
#' @examples
#' seqSet2Fasta(c("NC_031445.1", "NC_026892.1"), seq.label=c("accver", ", ", "orgname", ", chloroplast")) 
#'  
#' @rdname EUtils
seqSet2Fasta <- function(seqset.uid, by=5, folder=".", file.prefix="Refseq", file.extension="fasta", 
                         seq.label=c("accver", ", ", "orgname", ", chloroplast"), sleep=30) {
  fq <- c(seq(1, length(seqset.uid), by=by), (length(seqset.uid)+1))
  cat("break vector seqset.uid into", length(fq)-1, "subgroups, index = ", paste(fq, collapse = ", "), ".\n")
  
  require("reutils")
  f1 <- fq[1]
  for (f2 in fq[-1]) {
    uid <- as.character(seqset.uid[f1:(f2-1)])
    # efetch to get TSeqSet
    seqset <- efetch(uid, "nuccore", "fasta")
    seqset.df <- ComMA::parseTSeqSet(seqset$content)
    # create sequence label
    seqset.df$seq.name <- ">"
    for (sl in seq.label) {
      if (sl %in% colnames(seqset.df)) {
        seqset.df$seq.name <- paste0(seqset.df$seq.name, seqset.df[,sl])
      } else {
        seqset.df$seq.name <- paste0(seqset.df$seq.name, sl)
      }
    }

    fileConn <- file( file.path(folder, paste0(file.prefix, "-", f1, "-", (f2-1), ".", file.extension)) )
    cat("write", length(uid), "sequences", getFirstNInString(uid, 5), "to", fileConn, "\n")
    
    # covert 2 columns into a vector by rows
    seq.vec <- as.character(as.vector( t(as.matrix( seqset.df[,c("seq.name","sequence")] )) ))
    writeLines(seq.vec, fileConn)
    close(fileConn)
    
    Sys.sleep(sleep) # sleep 30 seconds
    f1 <- f2
  }
}
  

