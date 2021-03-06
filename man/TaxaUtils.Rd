% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TaxaUtils.R
\name{TaxaUtils}
\alias{TaxaUtils}
\alias{assignTaxaByRank}
\alias{combineTaxaAssign}
\alias{groupsTaxaMembers}
\alias{mergeCMTaxa}
\alias{prepTaxonomy}
\alias{subsetCM}
\alias{subsetTaxaTable}
\alias{summaryRank}
\alias{summaryTaxaAssign}
\title{Utils to preprocess taxa table}
\usage{
subsetTaxaTable(taxa.table, taxa.group = "assigned", rank = "kingdom",
  include = TRUE, ignore.case = TRUE)

subsetCM(community.matrix, taxa.table, taxa.group = NA, rank = NA,
  include = TRUE, ignore.case = TRUE, verbose = TRUE, drop.taxa = TRUE,
  merged.by = "row.names", ...)

prepTaxonomy(taxa.table, col.ranks = c("kingdom", "phylum", "class", "order",
  "family"), txt.unclassified = "unclassified", verbose = TRUE,
  pattern = "(\\\\s\\\\[|\\\\()(\\\\=|\\\\.|\\\\,|\\\\s|\\\\w|\\\\?)*(\\\\]|\\\\))")

mergeCMTaxa(community.matrix, taxa.table, classifier = c("MEGAN", "RDP"),
  min.conf = 0.8, has.total = 1, sort = TRUE, preprocess = TRUE,
  verbose = TRUE, mv.row.names = T,
  pattern = "(\\\\s\\\\[|\\\\()(\\\\=|\\\\.|\\\\,|\\\\s|\\\\w|\\\\?)*(\\\\]|\\\\))",
  col.ranks = c("kingdom", "phylum", "class", "order", "family"))

assignTaxaByRank(cm.taxa, unclassified = 0, aggre.FUN = sum,
  pattern = "(\\\\s\\\\[|\\\\()(\\\\=|\\\\.|\\\\,|\\\\s|\\\\w|\\\\?)*(\\\\]|\\\\))")

summaryTaxaAssign(ta.list, ta.OTU.list = list(), exclude.rank = c(-1),
  exclude.unclassified = TRUE, sort.rank = getRanks())

combineTaxaAssign(ta.list, keywords = c("Eukaryota"), ignore.case = TRUE,
  replace.to = c(), min.row.comb = 2)

summaryRank(ta.list, rank = "kingdom", exclude.unclassified = TRUE)

groupsTaxaMembers(taxa.assign, cm.taxa, rank = "phylum",
  rm.unclassified = TRUE, regex1 = "(\\\\|[0-9]+)", regex2 = "",
  ignore.case = TRUE, verbose = TRUE)
}
\arguments{
\item{taxa.table}{A data frame to contain taxonomic classifications of OTUs. 
Columns are taxonomy at the rank or lineage, rows are OTUs which need to 
match rows from community matrix. Use \code{\link{readTaxaTable}} to get it from file.}

\item{taxa.group}{The taxonomic group, the values can be 'all', 'assigned', or 
Group 'all' includes everything.
Group 'assigned' removes all uncertain classifications including 
'root', 'cellular organisms', 'No hits', 'Not assigned'. 
Alternatively, any high-ranking taxonomy in your taxonomy file 
can be used as a group or multi-groups (seperated by "|"), 
such as 'BACTERIA', 'Proteobacteria', etc. But they have to be 
in the same rank column in the file. Default to remove all 
uncertain classifications, even when group(s) assigned.}

\item{rank}{The rank to specify which column name in \code{taxa.table} to search.}

\item{include}{Define whether include or exclude given \code{taxa.group}. 
Default to TRUE.}

\item{ignore.case}{If TRUE, as default, case insensitive for taxon names.}

\item{community.matrix}{Community matrix (OTU table), where rows are 
OTUs or individual species and columns are sites or samples. See \code{\link{ComMA}}.}

\item{verbose}{More details. Default to TRUE.}

\item{drop.taxa}{TRUE, as default, to drop all taxonomy columns, 
and only keep \code{community.matrix} samples.}

\item{col.ranks}{A vector or string of column name(s) of taxonomic ranks in the taxa table, 
which will determine the aggregated abundence matrix. They have to be full set or subset of 
\code{c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")}. 
Default to \code{c("kingdom", "phylum", "class", "order", "family")}.}

\item{txt.unclassified}{The key word to represent unclassified taxonomy.}

\item{pattern}{The pattern for \code{\link{gsub}} "perl = TRUE". 
Default to "(\\s\\[|\\()(\\=|\\.|\\,|\\s|\\w|\\?)*(\\]|\\))".
Set NA to skip it.}

\item{classifier}{The classifier is used to generate \code{taxa.table}. 
Value is MEGAN or RDP. Default to MEGAN.}

\item{min.conf}{The confidence threshold to drop rows < \emph{min.conf}.}

\item{has.total}{If 0, then only return abundance by samples (columns) of community matrix. 
If 1, then only return total abundance. If 2, then return abundance by samples (columns) and total. 
Default to 1.}

\item{sort}{Sort the taxonomy rank by rank. Default to TRUE.}

\item{preprocess}{If TRUE, as default, replace 
"root|cellular organisms|No hits|Not assigned|unclassified sequences" from MEGAN result, 
or mark OTUs as 'unclassified' in RDP result whose confidence < \code{min.conf} threshold.}

\item{mv.row.names}{Default to TRUE to move the column 'Row.names' 
created by \code{\link{merge}} into data frame row.names, 
in order to keep the 1st column same as community matrix. 
Suggest not to change it.}

\item{cm.taxa}{The data frame combined community matrix with 
taxonomic classifications generated by \code{mergeCMTaxa}.
The row.names are OTUs, 1st column is the start of community matrix, 
\code{ncol.cm} column is the end of abundence, 
and \code{length(col.ranks)} columns taxonomy at different ranks. 
It should have attributes \code{ncol.cm} and \code{col.ranks}.

Note: From 1 to \code{ncol.cm} columns, the last column may be 'total' 
that is rowSums(cm) determined by \code{has.total} in \code{mergeCMTaxa}.}

\item{unclassified}{An interger to instruct how to deal with 
"unclassified" taxa. Default to 0, which keeps all "unclassified"
but moves them to the last rows. 
If 1, then remove the row whose taxon name is exact "unclassified". 
See the detail. 
If 2, then remove the row whose taxon name is exact "unclassified", 
but also merge all the rest "unclassified ???" to "unclassified rank",
such as "unclassified family".
If 3, then remove every rows containing "unclassified".
If 4, then do nothing.}

\item{aggre.FUN}{A function for \code{FUN} in \code{\link{aggregate}}. 
Default to \code{sum} to provide the reads abundance. 
Make \code{aggre.FUN=function(x) sum(x>0)} provide the OTU abundance.}

\item{ta.list, ta.OTU.list}{The list of taxonomic assignments 
created by \code{assignTaxaByRank} based on either number of reads or OTUs, 
where 'total' column is required. 
If ta.OTU.list is an empty list, as default, then do not count OTUs.
See \code{has.total} in \code{mergeCMTaxa} for the detail to get the total.}

\item{exclude.rank}{The first n elements (ranks) to exclude from the summary, default to -1, 
which is normally the kingdom.}

\item{exclude.unclassified}{Default to TRUE, not to count the taxonomy 
having the "unclassified" keyword.}

\item{sort.rank}{The order used to sort the summary dataframe by "rank" column, 
default to \code{c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")}.}

\item{keywords}{The vector of keywords for \code{\link{grep}}. 
The combined row will use the keyword as new row name.}

\item{replace.to}{The new names are used for combined rows, 
which should be either empty or the same length of the vector \code{keywords}.
If replace.to=c() as default, use the \code{keywords} as new row names.}

\item{min.row.comb}{The minimun number of rows from \code{\link{grep}} to combine. 
Default to 2, to ignore the single row selected by \code{\link{grep}} given a keyword.
Set to 1 to inlcude it in the combination process.}

\item{taxa.assign}{The data frame of taxonomic assignments with abundance
at the \code{rank}, where rownames are taxonomy at that rank, 
and columns are the sample names (may include total). It can be 
one element of the list generated by \code{assignTaxaByRank}. 
See the detail.}

\item{rm.unclassified}{Drop all unclassified rows (OTUs). Default to TRUE.}

\item{regex1, regex2}{Use for \code{\link{gsub}(regex1, regex2, row.names)} 
to remove or replace annotation from original labels. 
Default to \code{regex1="(\\|[0-9]+)", regex2=""}, 
which removes size annotation seperated by "|".}

\item{ignore.case}{Default to TRUE, same to \code{ignore.case} in \code{\link{grep}}.}

\item{rank}{The rank given to select the list of taxa assignments 
produced by \code{assignTaxaByRank}.}
}
\value{
\code{ncol.cm} and \code{col.ranks} are attributes of \code{cm.taxa} generated by \code{mergeCMTaxa}.

\code{ncol.cm} indicates how many column(s) is/are abundence in \code{cm.taxa}.

\code{col.ranks} records what ranks column(s) is/are in \code{cm.taxa}, 
which is also the input of \code{mergeCMTaxa}.
}
\description{
Utils to preprocess taxa table, 
and make it easy for visualization.
}
\details{
\code{subsetTaxaTable} takes or excludes a subset of given a taxa table at given rank.

\code{subsetCM} returns a subset community matrix 
regarding \code{taxa.group} at a given \code{rank} column 
in \code{taxa.table}, which is also the alternative choice 
of \code{mergeCMTaxa} if only simply \code{\link{merge}} is required. 
If either \code{taxa.group} or \code{rank} is NA, as default, 
then use the whole \code{taxa.table}, otherwise take the subset
of \code{taxa.table} by \code{subsetTaxaTable}.

\code{prepTaxonomy} replace repeated high rank taxa to 
unclassified high rank in MEGAN result, 
or replace the blank value to unclassified in RDP result,
in order to make taxonomy table \code{taxa.table} (can be \code{cm.taxa}) 
to make names look nice.
\code{col.ranks} vector have to be rank column names in taxa.table.

\code{mergeCMTaxa} creates a data frame \code{cm.taxa} combined community matrix with 
taxonomic classification table. The 1st column is "row.names" that are OTUs/individuals, 
the next "ncol.cm" columns are abundence that can be sample-based or total, 
and the last "length(col.ranks)" columns are the ranks. 

All sequences either classified as 
"root|cellular organisms|No hits|Not assigned|unclassified sequences"
from BLAST + MEGAN, or confidence < \emph{min.conf} threshold from RDP, 
are changed to "unclassified", which will be moved to the last row.

\code{assignTaxaByRank} provides a list of taxonomic assignments with abundance 
from community matrix at different rank levels, where rownames are taxonomy 
at that rank, and columns are the sample names (may include total). 
The function is iterated through \code{col.ranks}, and \code{\link{aggregate}}s 
abundance into taxonomy based on the rank in \code{col.ranks}.

\code{summaryTaxaAssign} summarises the number of reads, OTUs, and taxonomy
from the result of \code{assignTaxaByRank}.

\code{combineTaxaAssign} combines the total of taxonomy matching a given each of \code{keywords}
with the row names of the taxonomy assignment in the list from \code{assignTaxaByRank}.
The function is only working for the taxonomy assignment having 1 column "total" at the moment.

\code{summaryRank} directly converts the result of \code{assignTaxaByRank} 
at a given \code{rank} into a data frame as the summary.

\code{groupsTaxaMembers} groups the members (rows, also OTUs) from 
\code{cm.taxa} for each taxa in \code{taxa.assign} at the \code{rank}, 
and returns a list of members (OTUs) grouped by taxonomy. 
Default to drop all unclassified members (OTUs). 

It is impossible to trace back members after \code{assignTaxaByRank},
so that this function only has one option except the default, 
which assign the rest of members (OTUs) not picked up from other taxa  
into "unclassified". The result relies on using the identical \code{cm.taxa} 
in both \code{assignTaxaByRank} and \code{groupsTaxaMembers}.
}
\examples{
tt.sub <- subsetTaxaTable(tt.megan, taxa.group="Proteobacteria", rank="phylum")
tt.sub <- subsetTaxaTable(tt.megan, taxa.group="Cnidaria|Brachiopoda|Echinodermata|Porifera", rank="phylum", include=FALSE)

sub.cm <- subsetCM(cm, tt, taxa.group="BACTERIA", rank="kingdom")

tt <- prepTaxonomy(taxa.table, col.ranks=c("kingdom", "phylum", "class"))

cm.taxa <- mergeCMTaxa(community.matrix, tt.megan) 
ta.megan <- assignTaxaByRank(cm.taxa)

cm.taxa <- mergeCMTaxa(community.matrix, tt.rdp, classifier="RDP", has.total=0)
ta.rdp <- assignTaxaByRank(cm.taxa, unclassified=2)
colSums(ta.rdp[["phylum"]])

summary.ta.df <- summaryTaxaAssign(ta.list, ta.OTU.list)

combined.ta.list <- combineTaxaAssign(ta.list, c("Fungi", "Eukaryota", "Streptophyta|Viridiplantae", "Bacteria"))
combined.ta.list <- combineTaxaAssign(ta.list, c("Streptophyta|Viridiplantae"), replace.to=c("Plant"))

summary.kingdom.df <- summaryRank(ta.list, rank="kingdom")

taxa.members <- groupsTaxaMembers(ta.rdp[["phylum"]], tt.rdp)
taxa.members <- groupsTaxaMembers(ta.rdp[["family"]], tt.rdp, rank="family")

}
\keyword{taxonomy}

