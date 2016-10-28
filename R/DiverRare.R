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
        "replicates, using Jost diversities :", paste(levels, qs, sep="", collapse=", "), "\n.")
  
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
    print(summ)
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
  
  diver.rare.list <- list()
  for (ss in sample.sizes) {
    diver.rare.list[[paste("size.", ss, sep="")]] <- ComMA::getSubsamplesForDiverRare(cm, ss, replicates, levels=levels, qs=qs, 
                                                       progressBar=progressBar, verbose=verbose)
  }

  return(diver.rare.list)
}

getMultiDiverRare <- function(..., sample.sizes=10, replicates=10, levels=rep(c("gamma","alpha","beta"),3), 
                         qs=rep(0:2,each=3), progressBar=TRUE, min.sample.abundance=NA, 
                         non.rarefied.subsamples=TRUE, verbose=TRUE) {
                         
                         
}

plot.DiverRare <- function(diver.rare.list) {
  diver.rare.df <- data.frame(stringsAsFactors = F)
  
  for (i in 1:length(diver.rare.list)) {
    diver.rare <- diver.rare.list[[i]]
    diver.rare.df <- rbind(diver.rare.df, diver.rare$mean)
  }
  colnames(diver.rare.df) <- names(diver.rare.list[[1]]$mean)
  diver.rare.df$sample.size <- as.numeric(gsub("^.*?\\.","",names(diver.rare.list)))
  
  require(reshape2)
  melt.df <- melt(diver.rare.df, id="sample.size")
  melt.df[,"sample.size"] <- factor(melt.df[,"sample.size"], levels = sort(unique(melt.df[,"sample.size"])))
  
  gg.plot <- ComMA::ggLineWithPoints(melt.df[melt.df$variable=="beta1",], x.id="sample.size", y.id="value", group.id="variable", x.scale="discrete", 
                                     colour.id="variable", shape.id=shape.id, 
                                     point.data=point.data, line.or.point=line.or.point,
                                     line.type=line.type, line.alpha=line.alpha,
                                     title=title, x.lab=x.lab, y.lab=y.lab, 
                                     verbose=verbose, ...)
}

