# Author: Walter Xie
# Accessed on 22 Jun 2017

#' @name EUtils
#' @title EUtils NCBI
#'
#' @description Utils to process data downloaded from NCBI.
#' @import reutils,xml2
#' 
#' @details 
#' \code{parseTaxaSet} parses the result XML (DOCTYPE is TaxaSet) 
#' of \code{\link{efetch}} from taxonomy database into a data.frame, 
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
  as.data.frame(taxa.df)
}

#' @details 
#' \code{parseTSeqSet} parses the result XML (DOCTYPE is TSeqSet) 
#' of \code{\link{efetch}} from nuccore database into a data.frame, 
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
#' @export
#' @examples 
#' seqset <- efetch("NC_031445.1", "nuccore", "fasta")
#' seqset.df <- parseTSeqSet(seqset$content)
#' 
#' @rdname EUtils
parseTSeqSet <- function(xml) {
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
  }
  
  as.data.frame(seq.df)
}



