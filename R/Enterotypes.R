# enterotype  
# http://enterotype.embl.de/enterotypes_tutorial.sanger.R
#
# Modified by Walter Xie
# Accessed on 23 Nov 2016
#
#install.packages("cluster")
#install.packages("clusterSim")
#library(cluster)
#library(clusterSim)
#library(ade4)
#library(ggplot2)
#library(ellipse)


#' @name enterotype
#' @title Enterotype
#'
#' @description The enterotype method is introduced by
#' \url{http://enterotype.embl.de}. 
#' Besides the original functions copied from the website, 
#' we also add few more analyses to advancing this method, 
#' such as correlation between enterotypes and given/known groups.  
#' 
#' @details 
#' \code{enterotypePipeline} is a pipeline summarised from 
#' \url{http://enterotype.embl.de/enterotypes_tutorial.sanger.R}.
#' The steps of this pipeline are described below: 
#' 
#' 1. \code{getRelativeAbundance} turns taxa.assign into normalized probability distributions;
#' 
#' 2. \code{getJSD} calculates Jensen-Shannon divergence dissimilarity;
#' 
#' 3. \code{getClusters} lists Calinski-Harabasz (CH) Indices from 1 to n clusters;
#' 
#' 4. \code{data.cluster} picks up the optimised or selected k-cluster final result;
#' 
#' 5. \code{plotOptClusters} and \code{plotEnterotypes} make plots including BCA and PCoA.
#' 
#' ---------------------------------------------------------------------------------------
#' 
#' @keywords enterotype
#' @export
#' @examples 
#' enterotypePipeline(taxa.assign, n.max=20, fig.path="./figures", percent=0)
#'
#' @rdname enterotype
enterotypePipeline <- function(taxa.assign, n.max=20, k=NULL, fig.path=NA, percent=0.01, validate=TRUE, ...) {
  # 1. turn taxa.assign into normalized probability distributions
  relative.abund <- ComMA::getRelativeAbundance(taxa.assign)
  # 2. calculate Jensen-Shannon divergence dissimilarity
  jsd.dist <- ComMA::getJSD(relative.abund)
  # 3. Calinski-Harabasz (CH) Index of n clusters
  nclusters <- ComMA::getClusters(relative.abund, jsd.dist, n.max=n.max)
  # 4. the optimised or selected k-cluster final result
  data.cluster <- ComMA::getDataCluster(jsd.dist, k=k, nclusters=nclusters, validate=TRUE)
  
  # 5. plot
  ComMA::plotOptClusters(nclusters, fig.path=fig.path, n.max=n.max)
  ComMA::plotEnterotypes(relative.abund, jsd.dist, data.cluster, fig.path=fig.path, percent=percent, ...)
}

#' @details 
#' \code{getRelativeAbundance} returns a matrix of normalized probability distributions
#' of the abundance matrix for an input of \code{getJSD}, namely relative abundance.
#' 
#' @param taxa.assign The data frame of taxonomic assignments with abudence
#' at the \code{rank}, where rownames are taxonomy at that rank, 
#' and columns are the sample names. It can be 
#' one element of the list generated by \code{\link{assignTaxaByRank}}. 
#' @keywords enterotype
#' @export
#' @examples 
#' relative.abund <- getRelativeAbundance(taxa.assign)
#'
#' @rdname enterotype
getRelativeAbundance <- function(taxa.assign) {
  taxa.assign.prop <- prop.table(as.matrix(taxa.assign), 2)
}

#' @details 
#' \code{getJSD} return a distance \code{\link{dist}} object 
#' calculated by Jensen-Shannon divergence (JSD) metric.
#'  
#' @param data The matrix of normalized probability distributions 
#' of the abundance matrix, also called abundance distributions.
#' Rows are taxonomy at that rank and columns are samples. 
#' @keywords enterotype
#' @export
#' @examples 
#' jsd.dist <- getJSD(relative.abund)
#'
#' @rdname enterotype
getJSD <- function(data) {
  #nrow(data);ncol(relative.abund)
  #colSums(data)
  jsd.dist=dist.JSD(data)
  return(jsd.dist)
}

#' @details 
#' \code{getClusters} return the Calinski-Harabasz (CH) Index 
#' for choosig a number of clusters from 2 to \code{n.max}.
#'  
#' @param data.dist A distance \code{\link{dist}} object 
#' calculated by Jensen-Shannon divergence (JSD) metric.
#' @param n.max The number of clusters. Default to 20.
#' @keywords enterotype
#' @export
#' @examples 
#' nclusters <- getClusters(relative.abund, jsd.dist)
#'
#' @rdname enterotype
getClusters <- function(data, data.dist, n.max=20, ...) {
  if (n.max < 2)
    stop("The max number of clusters must >= 2!")
  
  require("clusterSim")
  #nclusters = index.G1(t(data), data.cluster, d = data.dist, centrotypes = "medoids")
  nclusters=NULL
  for (k in 1:n.max) { 
    if (k==1) {
      nclusters[k]=NA 
    } else {
      data.cluster_temp=pam.clustering(data.dist, k)
      nclusters[k]=index.G1(t(data), data.cluster_temp, d=data.dist, centrotypes = "medoids")
    }
  }
  return(nclusters)
}

#' @details 
#' \code{plotOptClusters} plots the Calinski-Harabasz (CH) Index 
#' for choosig a number of clusters from 2 to \code{n.max}.
#'  
#' @keywords enterotype
#' @export
#' @examples 
#' plotOptClusters(nclusters, fig.path="./figures", n.max=10)
#'
#' @rdname enterotype
plotOptClusters <- function(nclusters, fig.path=NA, n.max=20) {
  if (length(nclusters) < n.max)
    stop("n.max cannot excess length(nclusters) !")
  
  if (!is.na(fig.path))
    pdf(file.path(fig.path, "optimal-clusters-JSD.pdf"), width=6, height=6)
  
  # the optimal number of clusters
  plot(nclusters, type="h", xlab="number of clusters", ylab="CH index",
       main=paste0("Optimal number of clusters (Jensen-Shannon divergence)"), xaxt="n")
  axis(1, at = seq(1, n.max, by = 2))
  
  if (! is.na(fig.path))
    invisible(dev.off())   
}

#' @details 
#' \code{getDataCluster} returns the optimised or selected solution having 
#' \code{k} clusters assigned by \code{\link{pam}}.
#'  
#' @param k The number of clusters chosen for the optimised solution.
#' If NULL, the defult, it will take the number of clusters having 
#' the largest CH Index. 
#' @param nclusters The vector of nclusters, such as 
#' Calinski-Harabasz (CH) Index, to find the optimised solution.
#' @param validate Logical if it needs to validate. Default to TRUE.
#' @param as.vector Logical if TRUE, as default, to return a vector, 
#' otherwise return a \code{\link{pam.object}} representing the clustering. 
#' @keywords enterotype
#' @export
#' @examples 
#' data.cluster <- getDataCluster(jsd.dist, nclusters=nclusters)
#'
#' @rdname enterotype
getDataCluster <- function(data.dist, k=NULL, nclusters=NULL, as.vector=TRUE, validate=TRUE) {
  if (is.null(k)) {
    if (is.null(nclusters))
      stop("Please provide the nclusters, such as CH Index, to find the optimised solution !")
    k=match(max(nclusters[!is.na(nclusters)]), nclusters)
    cat("The optimised solution is ", k, " clusters.\n")
  }
  
  data.cluster=pam.clustering(data.dist, k=k, as.vector=as.vector)
  attr(data.cluster, "k") <- k
  if (as.vector) 
    return(data.cluster)
  
  #nclusters = index.G1(t(data), data.cluster, d = data.dist, centrotypes = "medoids")
  if (validate)
    validateCluster(data.cluster, data.dist)
  return(data.cluster)
}

#' @details 
#' \code{noise.removal} removes any row whose sum is samller than 
#' the given percentage of total sum of matrix. 
#' Advise to apply this function to data generated 
#' using short sequencing technologes, like Illumina or Solid.
#'  
#' @keywords enterotype
#' @export
#' @examples 
#' data=noise.removal(data, percent=0.1)
#'
#' @rdname enterotype
noise.removal <- function(data, percent=0.01){
  bigones <- rowSums(data)*100/sum(data) > percent 
  Matrix_1 <- data[bigones,]
  cat("percent = ", percent, ", remove", nrow(data)-nrow(Matrix_1), "rows.\n")
  return(Matrix_1)
}

#' @details 
#' \code{plotEnterotypes} is a mixed function to plot clusters, 
#' BCA \code{\link{bca}} and PCoA \code{\link{dudi.pco}}
#' from the optimised or selected solution.
#' 
#' Between-class analysis (BCA) was performed to support the 
#' clustering and identify the drivers for the enterotypes.
#' It is only available when \code{k} > 2.
#'  
#' @param data.cluster The clusters assigned by 
#' partitioning clustering \code{\link{pam}}.
#' @param attr.data A data frame to provide additional attributes 
#' for visualization. Default to give an empty data frame to do nothing.
#' @param percent The percentage threshold to remove the noise 
#' (ie. low abundant genera). Prior to the analysis. 
#' Default to 0.
#' @param addLabel Logical if add labels of points. Default to TRUE.
#' @param text.colour.id The column name in \code{attr.data} 
#' to colour texts of points.
#' @param fig.path The folder path to save figures. If NA, 
#' the default, do not plot the figure.
#' @keywords enterotype
#' @export
#' @examples 
#' plotEnterotypes(relative.abund, jsd.dist, data.cluster)
#'
#' @rdname enterotype
plotEnterotypes <- function(data, data.dist, data.cluster, attr.data=data.frame(), fig.path=NA, 
                            percent=0, addLabel=TRUE, text.colour.id=NULL, palette="Set1", postfix="",
                            plot.clusters=TRUE, plot.bca=TRUE, plot.pcoa=TRUE, 
                            cl.width=9, cl.height=6, width=8, height=8, verbose=TRUE, ...) {
  attr.data.cluster <- attributes(data.cluster)
  k <- as.numeric(attr.data.cluster$k)
  
  if (percent > 0)
    data=noise.removal(data, percent)
  
  if (plot.clusters)
    p.list <- ComMA::plotClusterAbundence(data, data.cluster, fig.path=fig.path, width=cl.width, height=cl.height)
  
  require(ade4)
  if (plot.bca && k > 2) {
    cat("Plot between-class analysis (BCA) for", k, "clusters.\n")
    
    #Between-class analysis (BCA) was performed to support the clustering and identify the drivers for the enterotypes. 
    obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=10)
    obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1) 
    cat("\nplot", k, "clusters between-class analysis (BCA).\n") 
    
    bca.df <- as.data.frame(obs.bet$ls)
    colnames(bca.df)[1:2] <- c("PC1", "PC2")
    bca.df$cluster <- data.cluster
    p.list[["BCA"]] <- plotSamples(bca.df, attr.data=attr.data, fig.path=fig.path, 
                                   addLabel=addLabel, text.colour.id=text.colour.id, palette=palette, 
                                   prefix=paste("enterotypes-bca", postfix, sep = "-"), 
                                   width=width, height=height, x.id="PC1", y.id="PC2", 
                                   title="Between Class Analysis", verbose=verbose, ...)
  }
  
  if (plot.pcoa) {
    cat("Plot principal coordinates analysis (PCoA) for", k, "clusters.\n")
    
    #principal coordinates analysis (PCoA) of a Euclidean distance matrix
    obs.pcoa=dudi.pco(data.dist, scannf=F, nf=3)
    cat("\nplot", k, "clusters principal coordinates analysis (PCoA).\n") 
    
    pcoa.df <- as.data.frame(obs.pcoa$li)
    colnames(pcoa.df)[1:3] <- c("PCo1", "PCo2", "PCo3")
    pcoa.df$cluster <- data.cluster
    p.list[["PCoA"]] <- plotSamples(pcoa.df, attr.data=attr.data, fig.path=fig.path, 
                                    addLabel=addLabel, text.colour.id=text.colour.id, palette=palette, 
                                    prefix=paste("enterotypes-pcoa", postfix, sep = "-"), 
                                    width=width, height=height, x.id="PCo1", y.id="PCo2", 
                                    title="Principal Coordinates Analysis", verbose=verbose, ...)
  }
  return(p.list)
}


#' @details 
#' \code{plotClusterAbundence} return a list of \code{\link{ggplot2}} objects
#' for relative abundance distribution in each cluster.
#'  
#' @keywords enterotype
#' @export
#' @examples 
#' p.list <- plotClusterAbundence(data, data.cluster)
#'
#' @rdname enterotype
plotClusterAbundence <- function(data, data.cluster, fig.path=NA, cluster.colours=c(), min.median=0, 
                                 x.lab="", y.lab="Relative abundence", width=9, height=6) {
  attr.data.cluster <- attributes(data.cluster)
  k <- attr.data.cluster$k
  
  if (!is(data, "matrix")) 
    stop("Data input has to be 'matrix' !")
  cat("\nPlot the relative abundance for", k, "clusters.\n")
  
  da.cl <- as.data.frame(cbind(t(data), cluster=data.cluster))
  require(reshape2)
  da.cl.melt <- melt(da.cl, id=c("cluster"))
  da.cl.melt <- da.cl.melt[da.cl.melt$value>0,]
  
  # TODO: cannot order by cluster and median in one data frame
  plot.list <- list()
  require(gg1L)
  for (cl in sort(unique(data.cluster))) {
    gg <- da.cl.melt[da.cl.melt$cluster==cl,]
    cat("Cluster", cl, "samples:", paste(rownames(da.cl[da.cl$cluster==cl,]), collapse = ", "), ".\n") 
    gg$cluster <- paste0("cluster", cl)
    
    #if (min.median > 0)
      
    gg$variable <- reorder(gg$variable, gg$value, median, order=TRUE)
    
    if (length(cluster.colours) > 0 && cl <= length(cluster.colours))
      cl.col <- cluster.colours[cl]
    else
      cl.col <- NULL
    
    p <- gg1L::ggBoxWhiskersPlot(gg, x.id="variable", y.id="value", fill.id="cluster",
                                  palette=cl.col, scale.type="manual",
                                  x.lab=x.lab, y.lab=y.lab, title=paste("Cluster", cl), 
                                  x.text.angle=90, no.legend="fill")
    plot.list[[cl]] <- p
    if(!is.na(fig.path))
      gg1L::pdf.ggplot(p, fig.path = file.path(fig.path, paste0("enterotype-",k,"-",cl,".pdf")), width=width, height=height)
  }
  return(plot.list)
}

######### correlation between enterotypes and known groups ##########
#' @details 
#' \code{corrEnterotypeToGroup} calculates Cramér's V between the enterotypes 
#' and the known groups (two categorical variables) from a same data set 
#' using Chi-Squared test \code{\link{chisq.test}}. 
#' \url{http://en.wikipedia.org/wiki/Cramér\%27s_V}.
#' 
#' @param group.id The column name in \code{attr.data} contains the known groups to 
#' compare with enterotypes.
#' @param simulate.p.value a logical indicating whether to compute p-values 
#' by Monte Carlo simulation, the default is FALSE.
#' @keywords enterotype
#' @export
#' @examples 
#' chi <- corrEnterotypeToGroup(jsd.dist, k=5, attr.data=env, group.id="land.use")
#'
#' @rdname enterotype
corrEnterotypeToGroup <- function(data.dist, k, attr.data, group.id, simulate.p.value=TRUE) {
  k.clusters <- pam.clustering(data.dist, k, as.vector=F)$clustering
  df.merge <- ComMA::mergeBy(as.data.frame(k.clusters), attr.data)
  
  clusters.id <- "k.clusters"
   # data frame bug: cannot use df.merge[, c(clusters.id, group.id)]
  pair.freq <- ComMA::freqUniqueValue(df.merge[, c(clusters.id, group.id)])
  pair.freq <- pair.freq[order(pair.freq[,clusters.id], pair.freq[,group.id]),]
  
  nrow=length(unique(pair.freq[,clusters.id])) 
  ncol=length(unique(pair.freq[,group.id]))
  tb.freq <- matrix(data=pair.freq$Freq, nrow=nrow, ncol=ncol, byrow=T)
  rownames(tb.freq) <- unique(pair.freq[,clusters.id])
  colnames(tb.freq) <- unique(pair.freq[,group.id])
  
  chi2v = chisq.test(tb.freq, simulate.p.value=simulate.p.value)
  #c(chi2v$statistic, chi2v$p.value)
  # V=sqrt( X-squared / (N * min(k-1, r-1)) ) 
  V=sqrt( as.numeric(chi2v$statistic) / ( sum(tb.freq) * min(ncol(tb.freq)-1, nrow(tb.freq)-1) ) )
  list(freq=tb.freq, chi=chi2v, V=V, p.value=chi2v$p.value)
}



######### internal ##########
# Jensen-Shannon divergence (JSD) (Endres & Schindelin, 2003)
# Kullback-Leibler divergence
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}
#data.dist=dist.JSD(data)

# Partitioning around medoids (PAM) clustering algorithm to cluster the abundance profiles. 
# PAM derives from the basic k-means algorithm, but has the advantage that it supports any 
# arbitrary distance measure and is more robust than k-means. It is a supervised procedure, 
# where the predetermined number of clusters is given as input to the procedure, which then 
# partitions the data into exactly that many clusters.
pam.clustering=function(x,k, as.vector=TRUE) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster=pam(as.dist(x), k, diss=TRUE)
  if (as.vector) {
    return(as.vector(cluster$clustering))
  }
  return(cluster)
}
#data.cluster=pam.clustering(data.dist, k=3)

# Observations with a large s(i) (almost 1) are very well clustered,  
# a small s(i) (around 0) means that the observation lies between two clusters,  
# and observations with a negative s(i) are probably placed in the wrong cluster.
validateCluster <- function(data.cluster, data.dist) {
  require(cluster)
  obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])
  cat("obs.silhouette =", obs.silhouette, "\n") 
}

# for both BCA and PCoA
plotSamples <- function(points.df, attr.data=data.frame(), fig.path=NA, addLabel=TRUE, 
                        colour.id="cluster", ellipsed.id="cluster", shape.id=NULL, text.colour.id=NULL, 
                        palette="Set1", prefix="enterotypes-pcoa-", width=8, height=8,
                        x.id="PCo1", y.id="PCo2", title="Principal Coordinates Analysis", verbose=TRUE, ...) {
  cluster.levels <- sort(unique(points.df$cluster))
  k <- length(cluster.levels)
  
  points.df$cluster <- factor(points.df$cluster, levels = cluster.levels)
  if (addLabel) {
    text.id="Row.names"
    if (nrow(attr.data) > 0) {
      points.df <- ComMA::mergeBy(points.df, attr.data, rm)
      #points.df[,] <- factor(points.df[,], levels = attr.levels)
    } else {
      points.df$Row.names <- rownames(points.df)
    }
  } else {
    text.id=NULL
  }
  require(gg1L)
  p <- gg1L::ggScatterPlot(points.df, x.id=x.id, y.id=y.id, colour.id=colour.id, ellipsed.id=ellipsed.id,
                            shape.id=shape.id, text.id=text.id, text.colour.id=text.colour.id, 
                            palette=palette, title=title, xintercept=0, yintercept=0, verbose=verbose, ...)
  gt <- gg1L::unclip.ggplot(p)  
  if(!is.na(fig.path))
    gg1L::pdf.gtable(gt, fig.path=file.path(fig.path, paste0(prefix, k,".pdf")), width=width, height=height) 
  return(p)
}

######### vegan cascadeKM ##########
plotCascadeKM <- function(fig.path, data, n.max) {
  library(vegan)
  
  ccas <- cascadeKM(t(data), 2, n.max)
  #ccas
  pdf(file.path(fig.path, paste0("optimal-clusters-cascadeKM.pdf")), width=6, height=6)
  plot(ccas, sortq=TRUE)
  invisible(dev.off()) 
  
  ccas <- cascadeKM(t(data), 2, 20, criterion = "ssi")
  pdf(file.path(fig.path, paste0("optimal-clusters-cascadeKM-ssi.pdf")), width=6, height=6)
  plot(ccas, sortq=TRUE)
  invisible(dev.off())
}

######### original plot fuction  ##########
plotEnterotypes.bak <- function(fig.path, data, data.dist, data.cluster, k, percent=0.01) {
  if (percent > 0)
    data=noise.removal(data, percent)
  
  require(ade4)
  ## plot 1
  obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=10)
  obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1) 
  
  pdf(file.path(fig.path, "enterotypes-between-class.pdf"), width=6, height=6)
  #  dev.new()
  s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F,sub="Between-class analysis")
  invisible(dev.off()) 
  
  #plot 2
  obs.pcoa=dudi.pco(data.dist, scannf=F, nf=k)
  
  pdf(file.path(fig.path, "enterotypes-principal-coordiante.pdf"), width=6, height=6)
  #  dev.new()
  s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F,sub="Principal coordiante analysis")
  invisible(dev.off())
}
