% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Enterotypes.R
\name{enterotype}
\alias{corrEnterotypeToGroup}
\alias{enterotype}
\alias{enterotypePipeline}
\alias{getClusters}
\alias{getDataCluster}
\alias{getJSD}
\alias{getRelativeAbundance}
\alias{noise.removal}
\alias{plotClusterAbundence}
\alias{plotEnterotypes}
\alias{plotOptClusters}
\title{Enterotype}
\usage{
enterotypePipeline(taxa.assign, n.max = 20, k = NULL, fig.path = NA,
  percent = 0.01, validate = TRUE, ...)

getRelativeAbundance(taxa.assign)

getJSD(data)

getClusters(data, data.dist, n.max = 20, ...)

plotOptClusters(nclusters, fig.path = NA, n.max = 20)

getDataCluster(data.dist, k = NULL, nclusters = NULL, as.vector = TRUE,
  validate = TRUE)

noise.removal(data, percent = 0.01)

plotEnterotypes(data, data.dist, data.cluster, attr.data = data.frame(),
  fig.path = NA, percent = 0, addLabel = TRUE, text.colour.id = NULL,
  palette = "Set1", postfix = "", plot.clusters = TRUE, plot.bca = TRUE,
  plot.pcoa = TRUE, cl.width = 9, cl.height = 6, width = 8,
  height = 8, verbose = TRUE, ...)

plotClusterAbundence(data, data.cluster, fig.path = NA,
  cluster.colours = c(), min.median = 0, x.lab = "",
  y.lab = "Relative abundence", width = 9, height = 6)

corrEnterotypeToGroup(data.dist, k, attr.data, group.id,
  simulate.p.value = TRUE)
}
\arguments{
\item{taxa.assign}{The data frame of taxonomic assignments with abudence
at the \code{rank}, where rownames are taxonomy at that rank, 
and columns are the sample names. It can be 
one element of the list generated by \code{\link{assignTaxaByRank}}.}

\item{n.max}{The number of clusters. Default to 20.}

\item{k}{The number of clusters chosen for the optimised solution.
If NULL, the defult, it will take the number of clusters having 
the largest CH Index.}

\item{fig.path}{The folder path to save figures. If NA, 
the default, do not plot the figure.}

\item{percent}{The percentage threshold to remove the noise 
(ie. low abundant genera). Prior to the analysis. 
Default to 0.}

\item{validate}{Logical if it needs to validate. Default to TRUE.}

\item{data}{The matrix of normalized probability distributions 
of the abundance matrix, also called abundance distributions.
Rows are taxonomy at that rank and columns are samples.}

\item{data.dist}{A distance \code{\link{dist}} object 
calculated by Jensen-Shannon divergence (JSD) metric.}

\item{nclusters}{The vector of nclusters, such as 
Calinski-Harabasz (CH) Index, to find the optimised solution.}

\item{as.vector}{Logical if TRUE, as default, to return a vector, 
otherwise return a \code{\link{pam.object}} representing the clustering.}

\item{data.cluster}{The clusters assigned by 
partitioning clustering \code{\link{pam}}.}

\item{attr.data}{A data frame to provide additional attributes 
for visualization. Default to give an empty data frame to do nothing.}

\item{addLabel}{Logical if add labels of points. Default to TRUE.}

\item{text.colour.id}{The column name in \code{attr.data} 
to colour texts of points.}

\item{group.id}{The column name in \code{attr.data} contains the known groups to 
compare with enterotypes.}

\item{simulate.p.value}{a logical indicating whether to compute p-values 
by Monte Carlo simulation, the default is FALSE.}
}
\description{
The enterotype method is introduced by
\url{http://enterotype.embl.de}. 
Besides the original functions copied from the website, 
we also add few more analyses to advancing this method, 
such as correlation between enterotypes and given/known groups.
}
\details{
\code{enterotypePipeline} is a pipeline summarised from 
\url{http://enterotype.embl.de/enterotypes_tutorial.sanger.R}.
The steps of this pipeline are described below: 

1. \code{getRelativeAbundance} turns taxa.assign into normalized probability distributions;

2. \code{getJSD} calculates Jensen-Shannon divergence dissimilarity;

3. \code{getClusters} lists Calinski-Harabasz (CH) Indices from 1 to n clusters;

4. \code{data.cluster} picks up the optimised or selected k-cluster final result;

5. \code{plotOptClusters} and \code{plotEnterotypes} make plots including BCA and PCoA.

---------------------------------------------------------------------------------------

\code{getRelativeAbundance} returns a matrix of normalized probability distributions
of the abundance matrix for an input of \code{getJSD}, namely relative abundance.

\code{getJSD} return a distance \code{\link{dist}} object 
calculated by Jensen-Shannon divergence (JSD) metric.

\code{getClusters} return the Calinski-Harabasz (CH) Index 
for choosig a number of clusters from 2 to \code{n.max}.

\code{plotOptClusters} plots the Calinski-Harabasz (CH) Index 
for choosig a number of clusters from 2 to \code{n.max}.

\code{getDataCluster} returns the optimised or selected solution having 
\code{k} clusters assigned by \code{\link{pam}}.

\code{noise.removal} removes any row whose sum is samller than 
the given percentage of total sum of matrix. 
Advise to apply this function to data generated 
using short sequencing technologes, like Illumina or Solid.

\code{plotEnterotypes} is a mixed function to plot clusters, 
BCA \code{\link{bca}} and PCoA \code{\link{dudi.pco}}
from the optimised or selected solution.

Between-class analysis (BCA) was performed to support the 
clustering and identify the drivers for the enterotypes.
It is only available when \code{k} > 2.

\code{plotClusterAbundence} return a list of \code{\link{ggplot2}} objects
for relative abundance distribution in each cluster.

\code{corrEnterotypeToGroup} calculates Cramér's V between the enterotypes 
and the known groups (two categorical variables) from a same data set 
using Chi-Squared test \code{\link{chisq.test}}. 
\url{http://en.wikipedia.org/wiki/Cramér\%27s_V}.
}
\examples{
enterotypePipeline(taxa.assign, n.max=20, fig.path="./figures", percent=0)

relative.abund <- getRelativeAbundance(taxa.assign)

jsd.dist <- getJSD(relative.abund)

nclusters <- getClusters(relative.abund, jsd.dist)

plotOptClusters(nclusters, fig.path="./figures", n.max=10)

data.cluster <- getDataCluster(jsd.dist, nclusters=nclusters)

data=noise.removal(data, percent=0.1)

plotEnterotypes(relative.abund, jsd.dist, data.cluster)

p.list <- plotClusterAbundence(data, data.cluster)

chi <- corrEnterotypeToGroup(jsd.dist, k=5, attr.data=env, group.id="land.use")

}
\keyword{enterotype}

