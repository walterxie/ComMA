# Plot prioritisation by ? diversity
# Author: Walter Xie
# Accessed on 3 Mar 2016

library(vegetarian)
library(reshape2)
library(ggplot2)

#' @title Plot prioritisation by Jost diversity
#'
#' @description Plot prioritisation by Jost diversity calculated from 
#' \pkg{vegetarian} \code{\link{d}}.
#' It uses a greedy algorithm to remove plots sequentially 
#' so as to minimize the loss of diversity among the remaining plots,
#' which always chooses the 1st plot if there are multi-results 
#' in each prioritisation loop.
#' Rank 1 is the most important plot and removed at the last, 
#' \emph{n} is the least important and removed in the beginning.
#' 
#' @param t.communityMatrix A transposed matrix from community matrix, 
#' where rows are plots (Use plots instead of subplots.), columns are OTUs.
#' @param lev Level of diversity to be calculated. Will accept: 'alpha,' 'beta,' or 'gamma.'
#' @param q Order of the diversity measure. Defaults to the Shannon case where q = 1.
#' @return 
#' A data frame with 2 columns: \emph{rank, diversity}. Rank 1 is the most important plot, 
#' \emph{n} is the least important, and row.names are plot names. For example,
#' \tabular{rrr}{
#'    \tab rank \tab diversity\cr
#'   CM30c39 \tab 28 \tab 1845.785714\cr
#'   CM30c44 \tab 27 \tab 1875.888889\cr
#'   CM31a5 \tab 26 \tab 1899.653846 
#' }
#' @export
#' @keywords plot prioritisation
#' @examples 
#' t.communityMatrix <- getCommunityMatrixT("16S", isPlot=TRUE, minAbund=1, taxa.group=""BACTERIA"")
#' priorPlots <- getPlotPriorByJostDiversity(t.communityMatrix, lev="gamma", q=1)
getPlotPriorByJostDiversity <- function(t.communityMatrix, lev, q){
	tmpCM <- t.communityMatrix
	removedSites <- c()

	jostDiversity <- d(t.communityMatrix, lev, q)
	print(paste("Original community Jost diversity ", lev, "(", q, ") = ", jostDiversity, " for ", nrow(t.communityMatrix), "sites", sep=""))
	# add orignal diversity
	maxRemainedDiversities <- c(jostDiversity)

	for (ra in nrow(t.communityMatrix):2) {
		if(ra != nrow(tmpCM)) stop(paste("incorrect nrow(tmpCM) ", nrow(tmpCM), " !=  ra ", ra, sep=""))
	
		maxRemainedDiv = 0
		removedSiteId = 0

		for (siteId in 1:nrow(tmpCM)) {
			tmpDiv <- d(tmpCM[-siteId,], lev, q)
			if (tmpDiv > maxRemainedDiv) {
				maxRemainedDiv = tmpDiv
				removedSiteId = siteId
			}
#			print(paste("siteId = ", siteId, ", remove siteId = ", removedSiteId, ", maxRemainedDiv = ", maxRemainedDiv, sep=""))  
		}

		if(maxRemainedDiv < 1 | removedSiteId < 1) stop(paste("incorrect max remained ", lev, "(", q, ") = ", maxRemainedDiv, ", when sites = ", nrow(tmpCM), sep=""))
	
		print(paste("remove site ", rownames(tmpCM)[removedSiteId], ", max remained ", lev, "(", q, ") = ", maxRemainedDiv, sep=""))  
		removedSites <- c(removedSites, rownames(tmpCM)[removedSiteId])
		maxRemainedDiversities <- c(maxRemainedDiversities, maxRemainedDiv)
	
		tmpCM <- tmpCM[-removedSiteId,]       
	}
	removedSites <- c(removedSites, rownames(tmpCM)[1]) # the last
	print(paste("Last site remained", rownames(tmpCM)[1]))
	
	if(length(removedSites) != length(maxRemainedDiversities)) 
		stop(paste("length of maxRemainedDiversities ", length(maxRemainedDiversities), " - 1 != sites ", length(maxDiv$rank), sep=""))

	# 1st removed (biggest rank) maxDiv$rank is the least important
	maxDiv <- data.frame(row.names=removedSites, rank=nrow(t.communityMatrix):1, diversity=maxRemainedDiversities) 
	maxDiv <- maxDiv[order(rownames(maxDiv)),]

	return(maxDiv)
}

#' Create a heat map for plot prioritisation result.
#' The plot prioritisation by Jost diversity can be 
#' calculated by \code{\link{getPlotPriorByJostDiversity}}. 
#' 
#' @param ranks.by.gourp A data frame contains ranks by groups for multi-dataset. 
#' Each column is the ranks of one dataset calculated by \code{\link{getPlotPriorByJostDiversity}}.
#' Each dataset can be grouped by genes, such as 16S, or by taxonomic groups, such as FUNGI. 
#' For example,
#' \tabular{rrrr}{
#'   plot \tab 16s \tab 18s \tab ITS\cr
#'   CM30c39 \tab 2 \tab 1 \tab 3\cr
#'   CM30c44 \tab 10 \tab 26 \tab 15\cr
#'   Plot01 \tab 6 \tab 5 \tab 6 
#' } 
#' @param fname.path The full path of image file.
#' @param title Graph title
#' @param y.lab The label of y-axis, such as plot names
#' @param low, high Refer to \pkg{ggplot2} \code{\link{scale_fill_gradient}}. Default to low="white", high="steelblue".
#' @param width, height Refer to \code{\link{pdf}}. Default to width=6, height=6.
#' @keywords plot prioritisation
#' @export
#' @examples 
#' ranks.by.gourp <- data.frame(plot=c("Plot03","Plot02","Plot01"), `16s`=c(3,2,1), `18s`=c(1,2,3), ITS=c(2,1,3), check.names = F)
#' ranks.by.gourp
#' heatmapPlotPrior(ranks.by.gourp, fname.path="plot-prior-example-heatmap.pdf")
heatmapPlotPrior <- function(ranks.by.gourp, fname.path, title="Plot Prioritisation Heatmap", y.lab="Plot", 
                             low="white", high="steelblue", width=6, height=6) {
  if (!is.element("plot", tolower(colnames(ranks.by.gourp))))
    stop("Data frame should have \"plot\" column !")
  
  breaks.rank <- round(seq(1, nrow(ranks.by.gourp), length.out = 5), digits = 0)
  ranks.melt <- melt(ranks.by.gourp, id=c("plot"))
  ranks.melt$plot <- factor(ranks.melt$plot, levels=unique(ranks.melt$plot))
  # variable is all group names, such as "16S" or "FUNGI"
  # value is ranks for each group
  p <- ggplot(ranks.melt, aes(x=variable, y=plot)) + geom_tile(aes(fill=value)) + 
    scale_fill_gradient(na.value="transparent", low=low, high=high, name="rank", breaks=breaks.rank) +
    ylab(y.lab) + ggtitle(title) +
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
  
  pdf(fname.path, width=width, height=height)	
  print(p)
  invisible(dev.off()) 
}


# maximum and minimum Jost diversity of all possible combinations of m_comb plots
# TODO old code, need to re-write 
getCombPlotPriorByJostDiversity <- function(t.communityMatrix, m_comb) { 
	plotsComb <- combn(rownames(t.communityMatrix), m_comb)
    d.comb <- matrix(0,nrow=9,ncol=ncol(plotsComb))
    rownames(d.comb) <- paste("$^",qs,"D_\\",levels,"$", sep="") #$^0D_\\gamma$
    
    for (i in 1:ncol(plotsComb)) {
        plots <- plotsComb[,i]
        matrixName <- paste(plots, collapse=" ")
    
        plotId <- which(rownames(t.communityMatrix) %in% plots)
        
        if (length(plotId) != m_comb) 
			stop(paste("Cannot find plot index from community matrix", plots))
    
		diversity_table <- diversityTable(t.communityMatrix[plotId,])

		for (j in  1:9) {
		   d.comb[j,i] <- unlist(diversity_table)[j]
		}
	}
	
	maxD <- matrix(0,nrow=9,ncol=(nrow(t.communityMatrix)+1))
	rownames(maxD) <- paste("$^",qs,"D_\\",levels,"$", sep="") #$^0D_\\gamma$
	colnames(maxD) <- c(rownames(t.communityMatrix), "Maximum")
	minD <- matrix(0,nrow=9,ncol=(nrow(t.communityMatrix)+1))
	rownames(minD) <- paste("$^",qs,"D_\\",levels,"$", sep="") #$^0D_\\gamma$
	colnames(minD) <- c(rownames(t.communityMatrix), "Minimum")
			
	for (j in  1:9) {
		maxD[j,ncol(maxD)] <- max(d.comb[j,])
		maxId <- which(d.comb[j,] == max(d.comb[j,]))
		plots <- plotsComb[,maxId]
		for (p in plots) {
			plotId <- which(colnames(maxD) %in% p)
			maxD[j,plotId] <- maxD[j,plotId] + 1
		}
		# convert frequency to ratio, last col is max value
		maxD[j,-ncol(maxD)] <- maxD[j,-ncol(maxD)]/max(maxD[j,-ncol(maxD)])  
		
		minD[j,ncol(minD)] <- min(d.comb[j,])
		minId <- which(d.comb[j,] == min(d.comb[j,]))
		plots <- plotsComb[,minId]
		for (p in plots) {
			plotId <- which(colnames(minD) %in% p)
			minD[j,plotId] <- minD[j,plotId] + 1
		}
		# convert frequency to ratio, last col is min value
		minD[j,-ncol(minD)] <- minD[j,-ncol(minD)]/max(minD[j,-ncol(minD)]) 
		
		print(paste("max", rownames(d.comb)[j], "index =", maxId, "; min index =", minId))
	}

#	maxD[,ncol(maxD)] <- formatC(signif(maxD[,ncol(maxD)],digits=5), digits=5,format="fg", flag="#")
#	minD[,ncol(minD)] <- formatC(signif(minD[,ncol(minD)],digits=5), digits=5,format="fg", flag="#")
#	maxD[maxD==0] <- ""
#	minD[minD==0] <- ""

	print(maxD)
	print(minD)

	list("maxD" = maxD, "minD" = minD)
}


# give same rank to the same remaining
# TODO in development
getPlotPriorByJostDiversity2 <- function(t.communityMatrix, lev, q){
  tmpCM <- t.communityMatrix[order(rownames(t.communityMatrix)),]
  
  ranks<-c()
  removedSites <- c()
  maxRemainedDiversities<-c()
  while(nrow(tmpCM) > 1) {		
    divs <- c()
    for (siteId in 1:nrow(tmpCM)) {
      tmpDiv <- d(tmpCM[-siteId,], lev, q)
      divs <- c(divs, tmpDiv)
      
      print(paste("siteId = ", siteId, ", tmpDiv = ", tmpDiv, sep="")) 
    }
    maxRemainedDiv = max(divs)	
    removedSiteId = which(divs==maxRemainedDiv)
    
    if(maxRemainedDiv < 0 | length(removedSiteId) < 1) 
      stop(paste("incorrect max remained ", lev, "(", q, ") = ", maxRemainedDiv, ", when sites = ", nrow(tmpCM), sep=""))
    
    print(paste("remove site ", rownames(tmpCM)[removedSiteId], ", max remained ", lev, "(", q, ") = ", maxRemainedDiv, sep=""))  
    
    ranks<-c(ranks, rep(nrow(tmpCM), length(removedSiteId)))
    removedSites <- c(removedSites, rownames(tmpCM)[removedSiteId])
    maxRemainedDiversities <- c(maxRemainedDiversities, rep(maxRemainedDiv, length(removedSiteId)))
    
    tmpCM <- tmpCM[-removedSiteId,]       
  }
  
  if (nrow(tmpCM) == 1) {
    ranks<-c(ranks, 1)
    removedSites <- c(removedSites, rownames(tmpCM)[1]) # the last
    maxRemainedDiversities <- c(maxRemainedDiversities, 0)
  }
  
  if(length(ranks) != length(maxRemainedDiversities)) 
    stop(paste("length of maxRemainedDiversities ", length(maxRemainedDiversities), " - 1 != sites ", length(ranks), sep=""))
  
  # 1st removed (biggest rank) maxDiv$rank is the least important
  maxDiv <- data.frame(row.names=removedSites, rank=ranks, diversity=maxRemainedDiversities)
  maxDiv <- maxDiv[order(rownames(maxDiv)),]
  
  return(maxDiv)
}
