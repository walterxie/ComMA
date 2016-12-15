# Plot prioritisation by ? diversity
# Author: Walter Xie
# Accessed on 3 Mar 2016

#' @name PlotPrioritisation
#' @title Plot prioritisation 
#'
#' @description Ranking of sample plots by their contributions 
#' to the total biodiversity.
#' 
#' @details \code{getPlotPrior.JostDiver} computes plot prioritisation 
#' by Jost diversity calculated from \pkg{vegetarian} \code{\link{d}}.
#' It uses a greedy algorithm to remove plots sequentially 
#' so as to minimize the loss of diversity among the remaining plots,
#' which always chooses the 1st plot if there are multi-results 
#' in each prioritisation loop.
#' Rank 1 is the most important plot and removed at the last, 
#' \emph{n} is the least important and removed in the beginning.
#' 
#' @param t.community.matrix A transposed matrix from community matrix, 
#' where rows are plots (Use plots instead of subplots.), columns are OTUs.
#' @param lev Level of diversity to be calculated. 
#' Will accept: 'alpha', 'beta', or 'gamma'.
#' @param q Order of the diversity measure. 
#' Defaults to the Shannon case where q = 1.
#' @param order.by How the result is ordered. 
#' Choose from 'sample', 'rank', or 'diversity'.
#' @return 
#' A data frame with 2 columns: \emph{rank, diversity}. 
#' Rank 1 is the most important plot, 
#' \emph{n} is the least important, and row.names are plot names. 
#' For example,
#' \tabular{rrr}{
#'    \tab rank \tab diversity\cr
#'   CM30c39 \tab 28 \tab 1845.785714\cr
#'   CM30c44 \tab 27 \tab 1875.888889\cr
#'   CM31a5 \tab 26 \tab 1899.653846 
#' }
#' @export
#' @keywords plot prioritisation
#' @examples 
#' plot.prior.g1 <- getPlotPrior.JostDiver(t.community.matrix, lev="gamma", q=1)
#' 
#' @rdname PlotPrioritisation
getPlotPrior.JostDiver <- function(t.community.matrix, lev=c("alpha","beta","gamma"), 
                                   q=1, order.by=c("sample","rank","diversity")){
  tmpCM <- t.community.matrix
  removedSites <- c()
  
  lev <- match.arg(lev)
  require(vegetarian)
  jost.diver <- d(t.community.matrix, lev=lev, q=q)
  cat("Original community", paste0(lev, "(", q, ")"), "=", jost.diver, 
      "for", nrow(t.community.matrix), "samples.\n")
  # add orignal diversity
  maxRemainedDiversities <- c(jost.diver)
  
  for (ra in nrow(t.community.matrix):2) {
    if(ra != nrow(tmpCM)) stop("incorrect nrow(tmpCM) ", nrow(tmpCM), " !=  ra ", ra, " !")
    
    maxRemainedDiv = 0
    removedSiteId = 0
    for (siteId in 1:nrow(tmpCM)) {
      tmpDiv <- d(tmpCM[-siteId,], lev=lev, q=q)
      if (tmpDiv > maxRemainedDiv) {
        maxRemainedDiv = tmpDiv
        removedSiteId = siteId
      }
      # print(paste("siteId = ", siteId, ", remove siteId = ", removedSiteId, ", maxRemainedDiv = ", maxRemainedDiv, sep=""))  
    }
    
    if(maxRemainedDiv < 1 | removedSiteId < 1) 
      stop("incorrect max remained ", lev, "(", q, ") = ", maxRemainedDiv, 
           ", when sites = ", nrow(tmpCM), " !")
    cat("remove sample", rownames(tmpCM)[removedSiteId], ": max remained", 
        paste0(lev, "(", q, ")"), "=", maxRemainedDiv, "\n")  
    
    removedSites <- c(removedSites, rownames(tmpCM)[removedSiteId])
    maxRemainedDiversities <- c(maxRemainedDiversities, maxRemainedDiv)
    # remove sample from cm
    tmpCM <- tmpCM[-removedSiteId,]       
  }
  removedSites <- c(removedSites, rownames(tmpCM)[1]) # the last
  cat("Last sample remained :", rownames(tmpCM)[1], "\n")
  
  if(length(removedSites) != length(maxRemainedDiversities)) 
    stop("length(maxRemainedDiversities) ", length(maxRemainedDiversities), 
         " - 1 != sites ", length(div.df$rank), " !")
  
  # 1st removed (biggest rank) div.df$rank is the least important
  div.df <- data.frame(row.names=removedSites, rank=nrow(t.community.matrix):1, 
                       diversity=maxRemainedDiversities) 
  div.df <- orderDiversityDF(div.df, order.by)
  return(div.df)
}

#' @details \code{getPlotPrior.PhyloAlpha} calculates plot prioritisation 
#' by phylogenetic alpha diversity from \code{\link{phylo.alpha}}.
#' It also can return the ranks based on species richness (SR), 
#' but they may be different to ranks calculated from 
#' \code{getPlotPrior.JostDiver} using gamma0 (also species richness).
#' 
#' @param phylo.tree,... The parameters passed to \code{\link{phylo.alpha}}.
#' @param taxa.match Logical, if taxa in phylogenies do not match OTUs in the community. 
#' If TRUE, as default, to use t.community.matrix and phylo.tree directly, 
#' otherwise to call \code{\link{match.phylo.comm}} 
#' before calculate phylogenetic alpha diversity.
#' @export
#' @keywords plot prioritisation
#' @examples 
#' phylo.alpha <- getPlotPrior.PhyloAlpha(t.community.matrix, phylo.tree)
#' 
#' @rdname PlotPrioritisation
getPlotPrior.PhyloAlpha <- function(t.community.matrix, phylo.tree, taxa.match=TRUE,
                                    order.by=c("sample","rank","diversity"), ...) {
  if (!taxa.match) {
    require(picante)
    combined <- match.phylo.comm(phylo.tree, t.community.matrix)
    pd.alpha <- ComMA::phylo.alpha(combined$comm, combined$phy, ...)
  } else {
    pd.alpha <- ComMA::phylo.alpha(t.community.matrix, phylo.tree, ...)
  }
  
  # add rank
  pd.alpha <- pd.alpha[order(pd.alpha$PD, decreasing = T),]
  pd.df <- data.frame(row.names=rownames(pd.alpha), rank=1:nrow(pd.alpha), diversity=pd.alpha$PD) 
  pd.df <- orderDiversityDF(pd.df, order.by)
  
  pd.alpha <- pd.alpha[order(pd.alpha$SR, decreasing = T),]
  sr.df <- data.frame(row.names=rownames(pd.alpha), rank=1:nrow(pd.alpha), diversity=pd.alpha$SR) 
  sr.df <- orderDiversityDF(sr.df, order.by)
  
  list(PD=pd.df, SR=sr.df)
}

#' @details \code{getPlotPrior} is a generic function including both 
#' \code{getPlotPrior.JostDiver} and \code{getPlotPrior.PhyloAlpha}, 
#' and it also handles multiple communities.
#' 
#' @param cm.list The list of community matrices. 
#' @param is.transposed If TRUE, then the community matrix is already
#' transposed to be the valid input of \code{\link{vegdist}}.  
#' Default to FASLE to transpose.
#' @param tre.list A list of phylo tree objects for 'pd.alpha' and 'sp.rich',
#' corresponding to \code{cm.list}. Default to an empty list. 
#' @param diversities The vector of diversities used to compute plot prioritisation.
#' The values are 'gamma0','gamma1','beta0','beta1','pd.alpha','sp.rich'. 
#' The first two are calculated by \code{\link{d}}, 
#' the last two by \code{\link{pd}}.
#' @export
#' @keywords plot prioritisation
#' @examples 
#' plot.prior.list <- getPlotPrior(cm.list, is.transposed=FALSE, diversities=c("gamma1","beta1"))
#' plot.prior.list <- getPlotPrior(cm.list, is.transposed=FALSE, tre.list=tre.list, diversities=c("gamma1","beta1","pd.alpha","sp.rich"))
#' @rdname PlotPrioritisation
getPlotPrior <- function(cm.list, is.transposed=FALSE, tre.list=list(), taxa.match=TRUE,
                         diversities=c("gamma0","gamma1","beta0","beta1","pd.alpha","sp.rich")) {
  cat("Plot prioritisation at", length(cm.list), "data sets", 
      ifelse(length(tre.list)>0, paste(", and", length(tre.list), "trees"), ""), ".\n") 
  # add names
  if (is.null(names(cm.list)))
    names(cm.list) <- 1:length(cm.list)
  
  plot.prior.list <- list()
  for (i in 1:length(cm.list)) {
    cm.name <- names(cm.list)[i]
    # if transposed
    if (!is.transposed)
      t.cm <- ComMA::transposeDF(cm.list[[i]])
    else
      t.cm <- cm.list[[i]]
    
    for (d in 1:length(diversities)) {
      div <- diversities[d]
      if (div=="gamma0") {
        plot.prior <- ComMA::getPlotPrior.JostDiver(t.cm, lev="gamma", q=0)
      } else if (div=="gamma1") {
        plot.prior <- ComMA::getPlotPrior.JostDiver(t.cm, lev="gamma", q=1)
      } else if (div=="beta0") {
        plot.prior <- ComMA::getPlotPrior.JostDiver(t.cm, lev="beta", q=0)
      } else if (div=="beta1") {
        plot.prior <- ComMA::getPlotPrior.JostDiver(t.cm, lev="beta", q=1)
      } else if (div=="pd.alpha" || div=="sp.rich") { 
        if(length(tre.list) < 1 || length(tre.list) != length(cm.list))
          stop("Invalid 'tre.list': a phylo tree object is required ", 
               "for each cm to caculate 'pd.alpha' or 'sp.rich' !")
        plot.prior.2 <- ComMA::getPlotPrior.PhyloAlpha(t.cm, tre.list[[i]], taxa.match=taxa.match)
        if (div=="pd.alpha")
          plot.prior <- plot.prior.2$PD
        else
          plot.prior <- plot.prior.2$SR
      } else {
        plot.prior <- NULL
        warning("Invalid diversity : ", div, " !")
      }
      plot.prior.list[[div]][[cm.name]] <- plot.prior
    }
  }
  plot.prior.list[["diversities"]] <- diversities
  plot.prior.list[["cm.names"]] <- names(cm.list)
  return(plot.prior.list)
}

#' @details \code{mergePlotPriorListOfDF} merges a list of data frames  
#' produced by \code{getPlotPrior} into one. Their columns must contain
#' 'rank' and 'diversity'.
#' 
#' @param plot.prior.list The output from \code{getPlotPrior}. 
#' @param suffixes The vector of suffixes added to distinguish colmun names 
#' merged by different data frames. Its length must equal to \code{plot.prior.list}.  
#' @export
#' @keywords plot prioritisation
#' @examples 
#' pp.df.list <- mergePlotPriorListOfDF(plot.prior.list)
#' 
#' @rdname PlotPrioritisation
mergePlotPriorListOfDF <- function(plot.prior.list, suffixes=c(), rm.prefix=TRUE) {
  diversities <- plot.prior.list[["diversities"]]
  if (length(suffixes) == 0)
    suffixes <- plot.prior.list[["cm.names"]]
  pp.df.list <- list()
  for (d in 1:length(diversities)) {
    div <- diversities[d]
    if (length(plot.prior.list[[div]]) != length(suffixes))
      stop("Invalid list : length(plot.prior.list[[div]]) != length(suffixes) !")
    pp.df <- ComMA::mergeListOfDF(plot.prior.list[[div]], suffixes=suffixes)
    rank.df <- pp.df[, grepl("^rank", colnames(pp.df))]
    div.df <- pp.df[, grepl("^diversity", colnames(pp.df))]
    if (rm.prefix) {
      colnames(rank.df) <- gsub("rank.", "", colnames(rank.df), ignore.case = T)
      colnames(div.df) <- gsub("diversity.", "", colnames(div.df), ignore.case = T)
    }
    pp.df.list[["rank"]][[div]] <- rank.df
    pp.df.list[["diversity"]][[div]] <- div.df
  }
  return(pp.df.list)
}

#' @details \code{plotPrioritisation} produces the list of \code{\link{ggHeatmap}}
#' given a list of output from \code{mergePlotPriorListOfDF}.
#' 
#' @param pp.df.list The output from \code{mergePlotPriorListOfDF}. 
#' @param at The values of heatmap, choose from "rank" or "diversity".  
#' @param row.id What row names are representing.
#' @export
#' @keywords plot prioritisation
#' @examples 
#' heatmap.list <- plotPrioritisation(pp.df.list, add.label=T)
#' 
#' @rdname PlotPrioritisation
plotPrioritisation <- function(pp.df.list, at=c("rank","diversity"), x.levels=c(), y.levels=c(), 
                               add.label=FALSE, label.digits=1, guide="colourbar", ...) {
  rank.df.list <-pp.df.list[[match.arg(at)]] 
  if (is.null(names(rank.df.list)))
    names(rank.df.list) <- 1:length(rank.df.list)
  
  heatmap.list <- list()
  for (i in 1:length(rank.df.list)) {
    pp.df <- rank.df.list[[i]]
    midpoint <- round(mean(c(max(pp.df), min(pp.df))))
    pp.df[,"samples"] <- rownames(pp.df)
    p.hm <- ComMA::ggHeatmap(pp.df, melt.id="samples", title="", x.levels=x.levels, y.levels=y.levels, 
                             add.label=add.label, label.digits=label.digits, guide=guide,
                             midpoint=midpoint, low="#f46d43", mid="#ffffbf", high="#3288bd" )
    heatmap.list[[names(rank.df.list)[i]]] <- p.hm
  }
  
  return(heatmap.list)
}

#' @details \code{plotPrioritisation.Phyloseq} produces 
#' one clusterd heatmap using \code{\link{phyloseq}} \code{\link{plot_heatmap}} 
#' and attaches an additional plot of a selected attribute  
#' from environmental meta-data, such as Elevation.
#' 
#' @param pp.df One element of the 'rank' list from \code{mergePlotPriorListOfDF}. 
#' @param attr.df Environmental meta-data. Rows are samples, 
#' and must include all rownames from \code{pp.df} 
#' @param y2.id,y2.lab The column of \code{attr.df} used 
#' to plot the 2nd figure on the right side, default to "Elevation".
#' @param grid.widths A unit vector giving the width of each two columns,
#' used by \code{\link{plot_grid}} in \pkg{cowplot} and default to \code{c(8, 2)}.
#' 
#' @note Use \code{theme_set(theme_bw(base_size=8))} to remove the grey backgroud.
#' @export
#' @keywords plot prioritisation
#' @examples 
#' hm.elv <- plotPrioritisation.Phyloseq(pp.df.list[["rank"]][[1]], env.plot, grid.widths = c(10,2))
#' plot(hm.elv$heatmap)
#' 
#' @rdname PlotPrioritisation
plotPrioritisation.Phyloseq <- function(pp.df, attr.df, y2.id="Elevation", y2.lab="Elevation (m)", 
                                         x.lab="Amplicon dataset", y.lab="Sample plot", 
                                        grid.widths = c(8, 2), ...) {
  if (! y2.id %in% colnames(attr.df) )
    stop("Invalid y2.id,", y2.id,  "not exsit in meta data column names !\n")
  
  midpoint <- mean(c(max(pp.df), min(pp.df)))
  require(phyloseq)
  phy <- phyloseq(otu_table(pp.df, taxa_are_rows = FALSE))
  ph <- plot_heatmap(phy, method = "RDA", distance = "euclidean") 
  require(ggplot2)
  ph <- ph + scale_fill_gradient2(low="#f46d43", mid="#ffffbf", high="#3288bd", midpoint = midpoint) + 
    xlab(x.lab) + ylab(y.lab) + coord_flip() +
    theme(legend.position="none", axis.text.x=element_text(size = 8, angle = -45, vjust = 1),
          axis.text.y = element_text(size = 8),
          panel.margin = unit(c(0,-0.5,0,0), "cm"))
  
  ### Extract data for secondary plot ###
  gb = ggplot_build(ph)
  p.ord <- as.data.frame(gb$data[[1]])
  # Sample col is defined by phyloseq
  p.ord$Sample <- gb$plot$data$Sample
  p.ord[,y2.id] <- attr.df[match(p.ord$Sample, rownames(attr.df)), y2.id]
  
  p2 <- ggplot(data = p.ord, aes_string("Sample", y2.id, group = "1")) + geom_line() + 
    coord_flip() + ylab(y2.lab) + xlab("") + 
    theme(axis.text.x = element_text(size = 8), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.margin = unit(c(0,0,0,-0.5), "cm"))
  
  require(cowplot)
  list(heatmap=plot_grid(ph, p2, align = 'h', rel_widths = grid.widths), 
       gg1=ph, gg2=p2)
}

#' @details \code{plotPrioritisation.NeatMap} produces 
#' one clusterd heatmap using \code{\link{NeatMap}} 
#' and attaches an additional plot of a selected attribute  
#' from environmental meta-data, such as Elevation.
#' 
#' @export
#' @keywords plot prioritisation
#' @examples 
#' hm.elv <- plotPrioritisation.NeatMap(pp.df.list[["rank"]][[1]], env.plot, grid.widths = c(10,2))
#' plot(hm.elv$heatmap)
#' 
#' @rdname PlotPrioritisation
plotPrioritisation.NeatMap <- function(pp.df, attr.df, y2.id="Elevation", y2.lab="Elevation (m)", 
                                       x.lab="Amplicon dataset", y.lab="Sample plot", 
                                       grid.widths = c(8, 3), ...) {
  if (! y2.id %in% colnames(attr.df) )
    stop("Invalid y2.id,", y2.id,  "not exsit in meta data column names !\n")
  
  midpoint <- mean(c(max(pp.df), min(pp.df)))
  require(ggplot2)
  require("NeatMap")
  p <- make.heatmap1(as.matrix(pp.df), 
                     column.method = "complete", column.metric="euclidean", # Column ordering method
                     row.method  ="complete", row.metric = "euclidean", 
                     #column.cluster.method = "none", # Column clustering method
                     #row.cluster.method = "none",
                     column.cluster.method = "complete", # Column clustering method
                     row.cluster.method = "complete",
                     column.labels = colnames(pp.df), row.labels = rownames(pp.df),
                     column.label.size = 2.5, row.label.size = 2.5) +
    scale_fill_gradient2(low="#f46d43", mid="#ffffbf", high="#3288bd", midpoint = midpoint) + 
    xlab(x.lab) + ylab(y.lab) + #coord_flip() +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(), legend.position="none",
          plot.margin = unit(c(0,2,2,0),"cm"), strip.text.x = element_text(colour = "red"))
  
  ### Get data for secondary plot ###
  gb = ggplot_build(p)
  #grid.draw(gb)
  
  #build$data[[4]]$y
  #build$data[[4]]$label
  #p.ord <- as.data.frame(gb$data[[2]]) # Without clusters
  p.ord <- as.data.frame(gb$data[[4]]) # With clusters
  #colnames(p.ord)[[1]] <- "Plot"
  p.ord <- p.ord[order(p.ord$y), ]
  p.ord$label <- factor(p.ord$label, levels = p.ord$label, ordered = TRUE)
  p.ord[,y2.id] <- attr.df[match(p.ord$label, rownames(attr.df)),y2.id]
  
  p2 <- ggplot(data = p.ord, aes_string("label", y2.id, group = "1")) + geom_line() + 
    coord_flip() + ylab(y2.lab) + xlab("") + 
    theme(axis.text.x = element_text(size = 8), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = unit(c(0,2,2,0),"cm"))
  
  require(cowplot)
  list(heatmap=plot_grid(p, p2, align = 'h', rel_widths = grid.widths), 
       gg1=p, gg2=p2)
}


############ internal #############
orderDiversityDF <- function(div.df, order.by=c("sample","rank","diversity")) {
  order.by <- match.arg(order.by)
  if (order.by=="rank") {
    div.df <- div.df[order(div.df$rank),]
  } else if (order.by=="diversity") {
    div.df <- div.df[order(div.df$rank, decreasing = T),]
  } else {
    div.df <- div.df[order(rownames(div.df)),]
  }
  return(div.df)
}


############ old code #############

# maximum and minimum Jost diversity of all possible combinations of m_comb plots
# TODO old code, need to re-write 
getCombPlotPrior.JostDiver <- function(t.community.matrix, m_comb) { 
  plotsComb <- combn(rownames(t.community.matrix), m_comb)
  d.comb <- matrix(0,nrow=9,ncol=ncol(plotsComb))
  rownames(d.comb) <- paste("$^",qs,"D_\\",levels,"$", sep="") #$^0D_\\gamma$
  
  for (i in 1:ncol(plotsComb)) {
    plots <- plotsComb[,i]
    matrixName <- paste(plots, collapse=" ")
    
    plotId <- which(rownames(t.community.matrix) %in% plots)
    
    if (length(plotId) != m_comb) 
      stop(paste("Cannot find plot index from community matrix", plots))
    
    diversity.table <- ComMA::diversityTable(t.community.matrix[plotId,])
    
    for (j in  1:9) {
      d.comb[j,i] <- unlist(diversity.table)[j]
    }
  }
  
  maxD <- matrix(0,nrow=9,ncol=(nrow(t.community.matrix)+1))
  rownames(maxD) <- paste("$^",qs,"D_\\",levels,"$", sep="") #$^0D_\\gamma$
  colnames(maxD) <- c(rownames(t.community.matrix), "Maximum")
  minD <- matrix(0,nrow=9,ncol=(nrow(t.community.matrix)+1))
  rownames(minD) <- paste("$^",qs,"D_\\",levels,"$", sep="") #$^0D_\\gamma$
  colnames(minD) <- c(rownames(t.community.matrix), "Minimum")
  
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
getPlotPrior.JostDiver2 <- function(t.community.matrix, lev, q){
  tmpCM <- t.community.matrix[order(rownames(t.community.matrix)),]
  
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
  
  # 1st removed (biggest rank) div.df$rank is the least important
  div.df <- data.frame(row.names=removedSites, rank=ranks, diversity=maxRemainedDiversities)
  div.df <- div.df[order(rownames(div.df)),]
  
  return(div.df)
}
