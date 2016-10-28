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
#' @export
#' @keywords rarefaction
#' @examples 
#' diver.rare <- getDiverRare(cm, replicates=5)
#' 
#' @rdname diverrare
getDiverRare <- function(cm, sample.sizes=10, replicates=10, levels=rep(c("gamma","alpha","beta"),3), 
                         qs=rep(0:2,each=3), progressBar=TRUE, min.sample.abundance=NA, verbose=TRUE) {
  require(ComMA)
  summary.cm <- ComMA::summaryCM.Vector(cm)
  rare.max <- summary.cm["max.sample.abun"]
  rare.min <- summary.cm["min.sample.abun"]
  
  if (!is.na(min.sample.abundance) && rare.min < min.sample.abundance) {
    cat("\nWarning: min sample abundance", rare.min, "< threshold", min.sample.abundance, 
        ", return NULL !\n")
    return(NULL)
  }
  if (length(sample.sizes) == 1) {
    points.mid <- round(sample.sizes/2)
    sample.sizes <- c(round(exp( seq(log(1), log(rare.min), length.out = points.mid)) )[1:(points.mid-1)], 
                      round(exp( seq(log(rare.min), log(rare.max), length.out = (sample.sizes-points.mid)) )) ) 
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


createRDPerSampleTable <- function(expId, isPlot, rmSingleton, taxa.group, pathFileStem, reps=10, min.size=100, verbose=T) {
  communityMatrix <- getCommunityMatrixT(expId, isPlot, rmSingleton, taxa.group)
  
  rare.max <- max(rowSums(communityMatrix))
  rare.min <- min(rowSums(communityMatrix))
  
  if (rare.min < min.size) {
    cat("\nWarning: min sample size", rare.min, "<", min.size, ", skip", taxa.group, "subset from", matrixNames[expId], ".\n")
    return(FALSE)
  }
  sample.sizes <- c(round(exp(seq(log(1), log(rare.min), length.out = 6)), digits = 0)[1:5], 
                    round(exp(seq(log(rare.min), log(rare.max), length.out = 6)), digits = 0))
  
  if (verbose) {
    cat("Subsampling: reps =", reps, ", min size allowed", min.size,".\n") 
    cat("Rarefaction :", matrixNames[expId], taxa.group, ", min =", rare.min, ", max =", rare.max, ".\n") 
  }
  
  alpha0MeanTable <- data.frame(row.names=rownames(communityMatrix), check.names=FALSE)
  alpha1MeanTable <- data.frame(row.names=rownames(communityMatrix), check.names=FALSE)
  
  for (ss in sample.sizes) {
    if (verbose)
      cat("sample size =", ss, ".\n") 
    
    rdMeanTable <- rrarefyPerSample(communityMatrix, ss, reps)
    alpha0MeanTable[,paste("size.", ss, sep="")] <- rdMeanTable$alpha0
    alpha1MeanTable[,paste("size.", ss, sep="")] <- rdMeanTable$alpha1
  }
  
  outputFile <- paste(pathFileStem, "rare", "alpha0", "table.csv", sep = "-") 
  write.csv(alpha0MeanTable, outputFile, row.names=TRUE, quote=FALSE)
  
  if (verbose) 
    cat("Write file", outputFile, ".\n")
  
  outputFile <- paste(pathFileStem, "rare", "alpha1", "table.csv", sep = "-")
  write.csv(alpha1MeanTable, outputFile, row.names=TRUE, quote=FALSE)
  
  if (verbose) 
    cat("Write file", outputFile, ".\n")
  
  return(TRUE)
}
