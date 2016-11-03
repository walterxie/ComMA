# Author: Walter Xie, Alexei Drummond
# Accessed on 28 Oct 2016
# Rarefaction using Jost diversities
 
#' @name diverrare
#' @title Rarefaction curves using Jost diversities
#'
#' @description Create rarefaction curves using Jost diversities.

#' @details \code{getSubsamplesForDiverRare} calculates subsamples for 
#' the rarefaction using Jost diversities, given a community matrix.
#' 
#' @param cm The community matrix.
#' @param is.transposed If TRUE, then \code{cm} is already
#' the transposed matrix of \code{\link{vegdist}} input. 
#' Default to FASLE.
#' @param sample.size An interger number representing subsample size for 
#' rarefying community using \code{\link{rrarefy}}.
#' It must be smaller than the minimun sample abundance, 
#' which can be calculated by \code{\link{summaryCM.Vector}}.
#' @param replicates The number of replicates to repeat subsampling.
#' @param levels The levels of Jost diversities in \code{\link{d}}.
#' Use \code{paste(levels, qs, sep="")} to get the actual measurements.
#' @param qs The order of Jost diversities in \code{\link{d}}.
#' @param progressBar Whether to print progress bar, default to TRUE.
#' @return 
#' \code{getSubsamplesForDiverRare} returns a list of results: 
#' \code{mean} is a named vector of means of each diversity measurement, 
#' \code{sterr} is the standard error of the mean, 
#' \code{subsamples} is a matrix of all subsampled points to calculate means and errors, 
#' where the columns are diversity measurements, rows are replicates.
#' @export
#' @keywords rarefaction
#' @examples 
#' diver.rare <- getSubsamplesForDiverRare(cm, 2000, 5)
#' 
#' @rdname diverrare
getSubsamplesForDiverRare <- function(cm, sample.size, replicates=10, levels=rep(c("gamma","alpha","beta"),3), 
                               qs=rep(0:2,each=3), is.transposed=FALSE, progressBar=TRUE, verbose=TRUE) {
  if (length(levels) != length(qs))
    stop("length(levels) != length(qs) !")
  
  require(ComMA)
  if (!is.transposed)
    cm <- ComMA::transposeDF(cm)
  
  size = length(levels)
  if (progressBar) {
    flush.console()
    pb <- txtProgressBar(min=1, max=replicates*size, style = 3)
  }
  if (verbose)
    cat("Subsampling rarefaction at", sample.size, "sample size and", replicates, 
        "replicates, using Jost diversities :", paste(levels, qs, sep="", collapse=", "), ".\n")
  
  rdt <- matrix(nrow=replicates,ncol=size)
  rownames(rdt) <- 1:replicates
  colnames(rdt) <- paste(levels, qs, sep="")
  require(vegan)
  for (i in 1:replicates) {
    cmrep <- rrarefy(cm, sample.size)
    for (j in 1:size) {
      if (progressBar) 
        setTxtProgressBar(pb, i*size + j)
      
      rdt[i,j] <- d(cmrep,lev=levels[j],q=qs[j])
    }
  }
  if (progressBar) 
    close(pb)
  
  mean.d <- apply(rdt,MARGIN=2, FUN=mean)
  sterr.d <- apply(rdt,MARGIN=2, FUN=sd)/sqrt(replicates)

  list(mean = mean.d, sterr = sterr.d, subsamples = rdt, sample.size = sample.size, 
       replicates = replicates, levels = levels, qs = qs)
}

#' @details \code{getDiverRare} returns a list of subsamples calculated 
#' from \code{getSubsamplesForDiverRare} given a community matrix.
#' Note: this community matrix must not be transposed.
#' 
#' @param sample.sizes The vector of subsample sizes for 
#' rarefying community using \code{\link{rrarefy}}.
#' It can be also an integer to define the length of vector, as default to 10, 
#' which allows to automatically create the vector of subsample sizes ranged by 
#' \code{min.sample.abundance} and \code{max.sample.abundance} of \code{cm}   
#' @param min.sample.abundance If \code{min.sample.abundance} of \code{cm} 
#' is smaller than this threshold, then return NULL.
#' @param non.rarefied.subsamples If TRUE, as default, the subsamples will include
#' non rarefied ones, whose \code{sample.size} > \code{min.sample.abundance}.
#' \code{\link{rrarefy}} will also create warnings.
#' @export
#' @keywords rarefaction
#' @examples 
#' diver.rare <- getDiverRare(cm, replicates=5)
#' 
#' @rdname diverrare
getDiverRare <- function(cm, sample.sizes=10, replicates=10, levels=rep(c("gamma","alpha","beta"),3), 
                         qs=rep(0:2,each=3), progressBar=TRUE, min.sample.abundance=NA, 
                         non.rarefied.subsamples=TRUE, verbose=TRUE) {
  require(ComMA)
  summary.cm <- ComMA::summaryCM.Vector(cm)
  rare.max <- summary.cm["max.sample.abun"]
  rare.min <- summary.cm["min.sample.abun"]
  if (verbose)
    print(summary.cm)
  if (!is.na(min.sample.abundance) && rare.min < min.sample.abundance) {
    cat("\nWarning: min sample abundance", rare.min, "< threshold", min.sample.abundance, 
        ", return NULL !\n")
    return(NULL)
  }
  if (length(sample.sizes) == 1) {
    if (non.rarefied.subsamples) {
      points.mid <- round(sample.sizes/2)
      sample.sizes <- c(round(exp( seq(log(1), log(rare.min), length.out = points.mid)) )[1:(points.mid-1)], 
                        round(exp( seq(log(rare.min), log(rare.max), length.out = (sample.sizes-points.mid+1)) )) ) 
    } else {
      sample.sizes <- round(exp( seq(log(1), log(rare.min), length.out = sample.sizes) ))
    }
  }
  if (length(sample.sizes) < 6)
    stop("Do not have enough points to plot rarefaction curve, please set bigger points (= ", 
         points, ") , or more sample.sizes (length = ", length(sample.sizes), ") !\n")
  
  diver.rare <- list()
  for (ss in sample.sizes) {
    diver.rare[[paste("size.", ss, sep="")]] <- ComMA::getSubsamplesForDiverRare(cm, ss, replicates, levels=levels, qs=qs, 
                                                       progressBar=progressBar, verbose=verbose)
  }
  return(diver.rare)
}

#' @details \code{getMultiDiverRare} returns a list of rarefaction curves  
#' calculated from \code{getDiverRare} given a list of community matrices.
#' Note: every community matrices must not be transposed.
#' 
#' @param cm.list The list of community matrices.
#' @param min.sample.abundance.list The list of \code{min.sample.abundance} 
#' for each \code{cm}, the default is NA.
#' @export
#' @keywords rarefaction
#' @examples 
#' multi.diver.rare <- getMultiDiverRare(cm.list, replicates=5)
#' 
#' @rdname diverrare
getMultiDiverRare <- function(cm.list, min.sample.abundance.list=NA, ...) {
  require(ComMA)
  multi.diver.rare <- list()
  names <- names(cm.list)
  for (i in 1:length(cm.list)) {
    cat("\nComputing rarefaction curves from community matrix", names[i], "...\n")
    
    if (is.na(min.sample.abundance.list)) {
      diver.rare <- ComMA::getDiverRare(cm.list[[i]], ...)
    } else {
      diver.rare <- ComMA::getDiverRare(cm.list[[i]], min.sample.abundance.list[[i]], ...)
    }
    if (is.null(names[i]) || is.na(names[i])) {
      multi.diver.rare[[i]] <- diver.rare
    } else {
      multi.diver.rare[[names[i]]] <- diver.rare
    }
  }                      
  return(multi.diver.rare)                       
}

#' @details \code{ggDiverRare} plots the rarefaction curves calculated 
#' from \code{getDiverRare} given a community matrix, 
#' or from \code{getMultiDiverRare} given a list of community matrices,
#' and returns a \code{\link{ggplot}} object.
#' 
#' @param diver.rare The rarefaction curves of one or more community matrices.
#' Use \code{multi.cm} to determine.
#' @param multi.cm If FALSE, as default, the \code{diver.rare} is calculated  
#' from one community matrix and no colour, 
#' otherwise it will be treated as a list of rarefaction curves
#' from a list of community matrices and coloured by \code{cm}.  
#' @param diversity The vector of diversity included in the graph, 
#' such as \code{c("gamma0","alpha0","beta1")}. The default is empty vector 
#' to inlcude all diversities in \code{diver.rare}. 
#' @export
#' @keywords rarefaction
#' @examples 
#' gg.plot <- ggDiverRare(diver.rare)
#' gg.plot <- ggDiverRare(multi.diver.rare, multi.cm=T, x.trans="log", auto.scale.x=T)
#' 
#' @rdname diverrare
ggDiverRare <- function(diver.rare, multi.cm=FALSE, diversity=c(),
                        line.or.point=3, point.size=2,
                        title="", x.lab="sample size", y.lab="diversity", ...) {
  diver.rare.df <- data.frame(stringsAsFactors = F)
  colour.id <- NULL
  group.id <- "variable"
  if (multi.cm) {
    colour.id <- "dataset"
    group.id <- colour.id
    diver.rare <- unlist(diver.rare, recursive=F)
  }
  
  for (i in 1:length(diver.rare)) {
    means <- diver.rare[[i]]$mean
    diver.rare.df <- rbind(diver.rare.df, means)
  }
  colnames(diver.rare.df) <- names(diver.rare[[1]]$mean)
  
  melt.id <- "sample.size"
  if (multi.cm) {
    melt.id <- c(melt.id, colour.id)
    # reomve size. prefix
    diver.rare.df[,melt.id[1]] <- as.numeric(gsub("^(.*)\\.size\\.(.*)","\\2",names(diver.rare)))
    # cm names
    diver.rare.df[,melt.id[2]] <- gsub("^(.*)\\.size\\.(.*)","\\1",names(diver.rare))
  } else {
    # reomve size. prefix
    diver.rare.df[,melt.id] <- as.numeric(gsub("^.*size\\.","",names(diver.rare)))
  }
  
  require(reshape2)
  melt.df <- melt(diver.rare.df, id=melt.id)
  # diversity filter to give different graph
  if (length(diversity) > 0)
    melt.df <- melt.df[melt.df$variable %in% diversity,]
  
  # assign x-axis
  melt.id <- "sample.size"
  melt.df[,melt.id] <- factor(melt.df[,melt.id], levels = sort(unique(melt.df[,melt.id])))
  # avoid Error: Discrete value supplied to continuous scale
  melt.df[,melt.id] =as.numeric(levels(melt.df[,melt.id]))[melt.df[,melt.id]]
#  if (multi.cm)
#    melt.df[,colour.id] <- factor(melt.df[,colour.id], levels = sort(unique(melt.df[,colour.id])))
  
  require(ComMA)
  gg.plot <- ComMA::ggLineWithPoints(melt.df, x.id=melt.id, y.id="value", 
                                     group.id=group.id, colour.id=colour.id, 
                                     point.size=point.size, line.or.point=line.or.point,
                                     x.scale="continuous", x.facet.id="variable", facet.scales="free",
                                     title=title, x.lab=x.lab, y.lab=y.lab, ...)
  return(gg.plot)
}

