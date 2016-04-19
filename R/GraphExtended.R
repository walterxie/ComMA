# Extended Graph 
# Author: Walter Xie
# Accessed on 15 Apr 2016

#' Use \code{\link{metaMDS}} in \code{\link{vegan}}
#' to create Nonmetric Multidimensional Scaling (NMDS) plot . 
#' 
#' @param comm Community data for NMDS plot. 
#' It is dissimilarities either as a \code{\link{dist}} structure 
#' or as a community data defined by \code{\link{vegan}} which is 
#' a transposed matrix of community matrix in \code{\link{ComMA}}.
#' Use \code{\link{transposeCM}} to rotate the data frame.
#' Detail to \code{\link{metaMDS}}.
#' @param distance Dissimilarity index used in \code{\link{vegdist}}.
#' @param attr.df A data frame to define how NMDS plot is coloured, shaped, or linked.
#' where its row names have to match row names of community data \code{comm}.
#' @keywords graph
#' @export
#' @examples
#' 
#' nmds.plot <- gtNMDSPlot(comm, attr.df, coloured.id="", shaped.id="", linked.id="")
gtNMDSPlot <- function(comm, attr.df, coloured.id, shaped.id, linked.id, 
                       distance = "bray", add.text=FALSE, ..., 
                       title="MDS", verbose=TRUE) {
  if (! missing(attr.df)) {
    if (! all(rownames(as.matrix(comm)) %in% rownames(attr.df)) )
      stop(paste("Invalid attr.df, rownames(as.matrix(comm)) should match rownames(attr.df)  !"))
  }

  agrs <- ""
  if (! missing(coloured.id)) {
    if (! coloured.id %in% colnames(attr.df))
      stop(paste("Invalid coloured.id,", coloured.id,  "not exsit in column names !"))
    agrs <- paste0(agrs, ", coloured.id=\"", coloured.id, "\"")
  } 
  if (! missing(linked.id)) { 
    if (! linked.id %in% colnames(attr.df) )
      stop(paste("Invalid linked.id,", linked.id,  "not exsit in column names !"))
    agrs <- paste0(agrs, ", linked.id=\"", linked.id, "\"")
  } 
  if (! missing(shaped.id)) { 
    if (! shaped.id %in% colnames(attr.df) )
      stop(paste("Invalid shaped.id,", shaped.id,  "not exsit in column names !"))
    agrs <- paste0(agrs, ", shaped.id=\"", shaped.id, "\"")
  }
    
  # Run metaMDS, get points and stress
  require(vegan)
  mds <- metaMDS(comm, distance = distance)
  pts.mds <- as.data.frame(mds$points)
  pts.mds <- pts.mds[order(rownames(pts.mds)),]
  
  if (title != "")
    title <- paste0(title, " (stress ", round(mds$stress, 2), ")")
  
  if (! missing(attr.df)) {
    rownames(pts.mds) <- tolower(rownames(pts.mds))
    rownames(attr.df) <- tolower(rownames(attr.df))
    
    pts.mds.merge <- merge(pts.mds, attr.df, by = "row.names")
    
    if (nrow(pts.mds.merge) != nrow(pts.mds) || nrow(pts.mds.merge) != nrow(attr.df)) 
      warning(paste("Some data points are missing after merge !", 
                    nrow(pts.mds.merge), "!=", nrow(pts.mds), "!=", nrow(attr.df) ))
  } else {
    pts.mds.merge <- pts.mds
    pts.mds.merge[,"Row.names"] <- rownames(pts.mds) 
  }
  
  if (add.text)
    agrs <- paste0(agrs, ", text.id=\"Row.names\"")
  
  if (verbose)
    cat("Passing additional argurments to ggScatterPlot : ", agrs, "\n")
  
  # Plot MDS ordination
  # ComMA::ggScatterPlot(pts.mds.merge, x.id="MDS1", y.id="MDS2", ..., title=title)
  ord.cmd = parse(text = paste0('ComMA::gtScatterPlot(pts.mds.merge, x.id="MDS1", y.id="MDS2",', 
                                agrs, ', ..., title=title)')) 
  g.d.gt <- eval(ord.cmd)
  
  return(g.d.gt)
}

#' Superimposed data ellipse on a scatter plot, which is extended from \code{\link{ggScatterPlot}}.
#' 
#' @param df.clusters The 3-column data frame for plotting. 
#' The first 2 columns are coordinates, and 3rd is cluster names. 
#' @param palette The palette to colour clusters using \code{\link{scale_colour_brewer}}. 
#' Default to 'Set1' (max 8 colours). Refer to \url{http://www.datavis.ca/sasmac/brewerpal.html}.
#' @return
#' A \code{\link{gtable}} object of scatter plot.
#' @keywords graph
#' @export
#' @examples 
#' df.clusters <- random2Clusters()
#' df.clusters
#' g.table <- gtScatterPlotEllipse(df.clusters, add.text=T)
#' pdfGtable(g.table, fig.path="clusters-scatter-plot.pdf") 
gtScatterPlotEllipse <- function(df.clusters, title="Clusters", palette="Set1", ...) {
  if (ncol(df.clusters) < 3)
    stop("Data frame should have 3 columns: first 2 columns are coordinates, 3rd is cluster names !")
  
  colnames(df.clusters)[1:3] <- c("PC1", "PC2", "cluster")
  # df.clusters$species <- paste(sapply(strsplit(rownames(df.clusters), "_"), "[[", 1), sapply(strsplit(rownames(df.clusters), "_"), "[[", 2), sep=".")
  df.clusters$cluster <- factor(df.clusters$cluster, levels = sort(unique(df.clusters$cluster)))

  g.d.gt <- ComMA::gtScatterPlot(df.clusters, x.id="PC1", y.id="PC2", coloured.id="cluster", 
                                 ellipsed.id="cluster", xintercept=0, yintercept=0, 
                                 palette=palette, title=title, ...)
  return(g.d.gt)
}


#' Percentage bar chart coloured by groups, which is extended from \code{\link{ggBarChart}}. 
#' 
#' @param df A data frame to \code{\link{melt}} and then make a percent bar chart. 
#' For example,
#' \tabular{rrrrr}{
#'   Phyla \tab 16s \tab 18s \tab ITS \tab TaxaGroup\cr
#'   Actinobacteria \tab 958 \tab 1 \tab 3 \tab Bacteria\cr
#'   Crenarchaeota \tab 1 \tab 0 \tab 0 \tab Archaea\cr
#'   Ascomycota \tab 2 \tab 765 \tab 971 \tab Fungi 
#' } 
#' @param melt.id A column name to \code{\link{melt}} and used to assign the colours.
#' @param title Graph title
#' @param x.lab,y.lab The label of x-axis or y-axis, such as plot names.
#' @param low, high Refer to \pkg{ggplot2} \code{\link{scale_fill_gradient}}. 
#' Default to low="white", high="steelblue".
#' @param autoWidth If TRUE, then use number of bars and legend columns 
#' to estimate pdf width automatically. Default to TRUE.
#' @keywords graph
#' @export
#' @examples 
#' taxa.phyla <- readFile("./data/examples/taxonomy97phyla.txt")
#' bar.chart <- ggPercentageBarChart(taxa.phyla, melt.id="TaxaGroup")
#' pdfGgplot(bar.chart$gg.plot, fig.path="taxa-percentage-bar.pdf", width=bar.chart$pdf.width, height=8) 
ggPercentageBarChart <- function(df, melt.id, title="Percent Bar Chart", x.lab="", y.lab="", autoWidth=TRUE) {
  if (!is.element(tolower(melt.id), tolower(colnames(df))))
    stop(paste0("Data frame column names do NOT have \"", melt.id, "\" for melt function !"))
  
  require(reshape2)
  df.melt <- melt(df, id=c(melt.id))
  #df.melt[,"variable"] <- factor(df.melt[,"variable"], levels = sort(unique(df.melt[,"variable"])))
  
  # move unclassified group to the last of legend 
  legend.ord <- as.character(unique(df[,melt.id]))
  id.match <- grep("unclassified", legend.ord, ignore.case = TRUE)
  if (length(id.match) > 0)
    legend.ord <- legend.ord[c(setdiff(1:length(legend.ord), id.match),id.match)]
  df.melt[,melt.id] <- factor(df.melt[,melt.id], levels = rev(legend.ord))
  
  pale <- ComMA::getMyPalette(length(legend.ord))
  if (length(legend.ord) > length(pale)) {
    require(colorspace)
    pale <- rainbow_hcl(length(legend.ord))
  }
  # number of columns for legend
  legend.col = ceiling(length(legend.ord) / 25)
  
  p <- ComMA::ggBarChart(df.melt, x.id="variable", y.id="value", fill.id=melt.id, bar.pos="fill", 
                  y.scale="per", title=title, x.lab=x.lab, y.lab=y.lab, palette=pale, legend.col=legend.col)
  
  if (autoWidth)
    pdf.width = 1 + legend.col*2.5 + length(unique(df.melt[,"variable"])) * 0.2
  
  # Return a list containing the filename
  list(
    pdf.width = pdf.width,
    legend.len = length(legend.ord), # the length of legend
    legend.ord = legend.ord, # the legend in the order
    gg.plot = p # ggplot
  )
}
