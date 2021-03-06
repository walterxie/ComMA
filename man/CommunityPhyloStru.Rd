% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CommunityPhyloStrut.R
\name{CommunityPhyloStru}
\alias{CommunityPhyloStru}
\alias{phylo.alpha}
\alias{phylo.beta.dist}
\alias{phylo.mntd}
\alias{phylo.mpd}
\alias{phylo.unifrac.dist}
\alias{plotPDBeta}
\alias{printResult}
\title{Community Phylogenetic Structure from \pkg{picante} Package}
\usage{
phylo.alpha(t.community.matrix, phylo.tree, include.root = TRUE,
  ORD.RES = function(res) {     res[order(rownames(res)), ] },
  verbose = TRUE)

phylo.mpd(t.community.matrix, phylo.tree, ORD.RES = function(res) {    
  res[order(rownames(res)), ] }, verbose = TRUE)

phylo.mntd(t.community.matrix, phylo.tree, ORD.RES = function(res) {    
  res[order(rownames(res)), ] }, verbose = TRUE)

phylo.beta.dist(t.community.matrix, phylo.tree, verbose = TRUE)

phylo.unifrac.dist(t.community.matrix, phylo.tree, verbose = TRUE)

printResult(pd.alpha = NULL, pd.mpd = NULL, pd.mntd = NULL,
  pd.beta.dist = NULL, filePath, fileStem, tableFile = NULL)

plotPDBeta(pd.beta.dist, title = "Phylogenetic beta diversity", xlab = "",
  sub = "", ...)
}
\arguments{
\item{t.community.matrix}{A transposed matrix from community matrix, 
where rows are samples, columns are OTUs.}

\item{phylo.tree}{A rooted tree of phylo object}

\item{ORD.RES}{The function how to order the sample names in the result}

\item{verbose}{default TRUE}

\item{pd.alpha}{The result of phylogenetic alpha diversity, ignore it if NULL.}

\item{pd.mpd}{The result of MPD, ignore it if NULL.}

\item{pd.mntd}{The result of MNTD, ignore it if NULL.}

\item{pd.beta.dist}{The result of phylogenetic beta diversity, ignore it if NULL.}

\item{filePath}{The path to create the file if tableFile is not NULL.}

\item{fileStem}{The file prefix to identify the result multi-conditions from multi-datasets}

\item{tableFile}{If NULL, then print the results to console, otherwise print them to the file. Default to NULL.}
}
\description{
A wrapper class using \pkg{picante} Package. 
Data input \strong{t.community.matrix} is 
a transposed matrix from community matrix we defined in \pkg{ComMA}.
\strong{phylo.tree} is a rooted tree of phylo object, 
which can get from \pkg{ape} \code{\link{read.tree}}.
}
\details{
\code{phylo.alpha} returns a data frame of the PD and species richness (SR) 
values for all samples from \code{\link{pd}} in \pkg{picante}. 
Phylogenetic alpha diversity (PD) index is proposed by Faith (1992).

\code{phylo.mpd} returns \pkg{picante} \code{\link{ses.mpd}}, 
which is MPD standardized effect size of mean pairwise distances in communities.
When used with a phylogenetic distance matrix, equivalent to -1 times the Nearest Taxon Index (NTI).

\code{phylo.mntd} returns \pkg{picante} \code{\link{ses.mntd}}, 
which is MNTD standardized effect size of mean nearest taxon distances in communities.

\code{phylo.beta.dist} returns \code{\link{dist}} object from \pkg{picante} \code{\link{comdist}}, 
which is phylogenetic beta diversity (Steven Kembel).

\code{phylo.unifrac.dist} calculates unweighted UniFrac and 
returns dist object from \pkg{picante} \code{\link{unifrac}}.

\code{printResult} prints intermediate data to either file or console.

\code{plotPDBeta} plots the dendrogram of phylogenetic beta diversity (Steven Kembel)
}
\note{
If taxa in phylogenies do not match OTUs in the community, 
then use \code{\link{match.phylo.comm}} before \code{phylo.alpha}. 
See examples. But prefer \code{comm > phy} in practical, 
otherwise \code{combined$phy} will be unrooted tree. 
Therefore, in the case of single-species samples the PD will 
be equal to NA (include.root=FALSE), see \code{\link{pd}}.
}
\examples{
pd.alpha <- phylo.alpha(t.community.matrix, phylo.tree)
pd.alpha
\\tabular{rrr}{
   \\tab PD \\tab SR\\cr
  CM30b51 \\tab 402.0822 \\tab 2796\\cr
  CM30b58 \\tab 456.4397 \\tab 2634\\cr
  CM30b60 \\tab 554.9873 \\tab 2931
}

combined <- match.phylo.comm(phylo.tree, t.communityMatrix)
pd.alpha <- phylo.alpha(combined$comm, combined$phy, ...) 

pd.mpd <- phylo.mpd(t.community.matrix, phylo.tree)

pd.mntd <- phylo.mntd(t.community.matrix, phylo.tree)
pd.beta.dist <- phylo.beta.dist(t.community.matrix, phylo.tree)

unifrac.dist <- phylo.unifrac.dist(t.community.matrix, phylo.tree)

tableFile <- file.path(workingPath, "report.tex")
pdFilePath <- file.path(workingPath, "data", "pd")
mkdir(pdFilePath) 
fileStem <- paste("16S", postfix(taxa.group, isPlot, sep="-"), sep = "-")
printResult(pd.alpha, pd.mpd, pd.mntd, pd.beta.dist, pdFilePath, fileStem, tableFile)

library(pryr);library(gg1L)
p \%<a-\% plotPDBeta(pd.beta.dist)
gg1L::pdf.plot(p, fig.path="pd-beta.pdf", width=10, height=5)  

}

