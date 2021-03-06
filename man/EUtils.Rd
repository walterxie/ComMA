% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EUtils.R
\name{EUtils}
\alias{EUtils}
\alias{batchDownload}
\alias{parseTSeqSet}
\alias{parseTaxaSet}
\alias{seqSet2Fasta}
\title{EUtils NCBI}
\usage{
parseTaxaSet(xml)

parseTSeqSet(xml, save.seq = T)

batchDownload(uid.vec, db = NULL, rettype = NULL, retmode = NULL,
  out.file = "res.txt", sleep = 10, ...)

seqSet2Fasta(seqset.uid, by = 5, folder = ".", file.prefix = "Refseq",
  file.extension = "fasta", seq.label = c("accver", ", ", "orgname",
  ", chloroplast"), sleep = 30)
}
\arguments{
\item{xml}{The XML result of \code{\link{efetch}}.}

\item{save.seq}{If TRUE as default, save sequences 
into data.frame returned by \code{parseTSeqSet}.}

\item{uid.vec}{The vector of \code{uid} for \code{\link{efetch}} to NCBI data. 
The \code{uid} will either look like NC_016668.1 ("accver") or KM462867 ("INSDC").}

\item{db, rettype, retmode, ...}{The arguments of \code{\link{efetch}}.
See \url{https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly}
for all supported databases and their available values.}

\item{out.file}{The file to write all results directly without parsing.
Note: if \code{length(uid) > 500}, \code{outfile} of \code{\link{efetch}} is required.}

\item{sleep}{Please be nice to give enough time to break. Default 10 seconds.}

\item{seqset.uid}{The vector of \code{uid} for \code{\link{efetch}} to download sequence,
where db = "nuccore", rettype = "fasta". 
The \code{uid} will either look like NC_016668.1 ("accver") or KM462867 ("INSDC").}

\item{by}{Split \code{seqset.uid} into subgroups by the given number, defaul to 20.}

\item{folder, file.prefix, file.extension}{Determine the file name.}

\item{seq.label}{The vector of string to determine how to label the sequence. 
If the element is one of the column name data.frame from \code{parseTSeqSet}, 
then the label of that position will be the value of that column. 
The available columns are "seqtype", "gi", "accver", "taxid", "orgname", "defline", "length".}
}
\description{
Utils to process data downloaded from NCBI.
}
\details{
\code{parseTaxaSet} parses the taxonomy XML (DOCTYPE is TaxaSet) 
as the result of \code{\link{efetch}} from taxonomy database into a data.frame, 
which inlcudes "TaxId", "ScientificName", "Rank", "Lineage", "Division", 
and the format of taxa.table from "kingdom" to "genus".

<TaxaSet>
<Taxon>
  <TaxId>123685</TaxId>
  <ScientificName>Oryzias minutillus</ScientificName>
  <ParentTaxId>8089</ParentTaxId>
  <Rank>species</Rank>
  <Division>Vertebrates</Division>
  <Lineage>cellular organisms; Eukaryota; Opisthokonta; Metazoa; ...</Lineage>
  <LineageEx>
   <Taxon>
    <TaxId>131567</TaxId>
    <ScientificName>cellular organisms</ScientificName>
    <Rank>no rank</Rank>
   </Taxon>
   ...  
  </LineageEx>
</Taxon>
...
</TaxaSet>

\code{parseTSeqSet} parses the TinySeq XML (DOCTYPE is TSeqSet) 
as the result of \code{\link{efetch}} from nuccore database into a data.frame, 
which inlcudes "TaxId", "ScientificName", "ACCESSION", "Lineage", "sequence".

<TSeqSet>
  <TSeq>
    <TSeq_seqtype value="nucleotide"/>
    <TSeq_gi>1079489517</TSeq_gi>
    <TSeq_accver>NC_031445.1</TSeq_accver>
    <TSeq_sid>gnl|NCBI_GENOMES|60824</TSeq_sid>
    <TSeq_taxid>126358</TSeq_taxid>
    <TSeq_orgname>Abeliophyllum distichum</TSeq_orgname>
    <TSeq_defline>Abeliophyllum distichum chloroplast, complete genome</TSeq_defline>
    <TSeq_length>155982</TSeq_length>
    <TSeq_sequence>CATTTTAGTTATGGGC...GCTGT</TSeq_sequence>
  </TSeq>
</TSeqSet>

\code{batchDownload} downloads NCBI data given a vector of \code{uid} 
using \code{\link{efetch}} one at a time for each \code{uid}, 
and writes all results directly to a file without parsing.

\code{seqSet2Fasta} downloads reference sequences given their \code{uid} 
using \code{\link{efetch}}. \code{parseTSeqSet} is used to parse the 
result of \code{\link{efetch}}.
}
\examples{
library("reutils")
taxa <- efetch(c("123685", "8089", "8088"), "taxonomy")
taxa.df <- parseTaxaSet(taxa$content)

seqset <- efetch("NC_031445.1", "nuccore", "fasta")
seqset.df <- parseTSeqSet(seqset$content)

batchDownload(c("NC_031445.1", "NC_026892.1"), "nuccore", "gb", out.file="res.gbff") 
 
seqSet2Fasta(c("NC_031445.1", "NC_026892.1"), seq.label=c("accver", ", ", "orgname", ", chloroplast")) 
 
}
\keyword{eutils}

