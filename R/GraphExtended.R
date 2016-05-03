# Extended Graph 
# Author: Walter Xie
# Accessed on 15 Apr 2016

#' NMDS Plot
#' 
#' Use \code{\link{metaMDS}} in \code{\link{vegan}}
#' to create Nonmetric Multidimensional Scaling (NMDS) plot.
#' It is extended from \code{\link{gtScatterPlot}}. 
#' 
#' @param comm Community data for NMDS plot. 
#' It is dissimilarities either as a \code{\link{dist}} structure 
#' or as a community data defined by \code{\link{vegan}} which is 
#' a transposed matrix of community matrix in \code{\link{ComMA}}.
#' Use \code{\link{transpose.df}} to rotate the data frame.
#' Detail to \code{\link{metaMDS}}.
#' @param distance Dissimilarity index used in \code{\link{vegdist}}.
#' @param attr.df A data frame to define how NMDS plot is coloured, shaped, or linked.
#' where its row names have to match row names of community data \code{comm}.
#' @keywords graph
#' @export
#' @examples
#' 
#' nmds.plot <- gtNMDSPlot(comm, env, colour.id="FishSpecies", shape.id="FeedingPattern", add.text=T)
#' require(grid)
#' grid.draw(nmds.plot)
gtNMDSPlot <- function(comm, attr.df, colour.id=NULL, shape.id=NULL, link.id=NULL, 
                       add.text=TRUE, text.data = NULL, text.size=3, palette=NULL,
                       xintercept=NULL, yintercept=NULL, 
                       distance="bray", title="MDS", verbose=TRUE, ...) {
  if (! missing(attr.df)) {
    if (! all(rownames(as.matrix(comm)) %in% rownames(attr.df)) )
      stop("Invalid attr.df,", paste(rownames(as.matrix(comm)), collapse = ","), 
           "should match", paste(rownames(attr.df), collapse = ","), "!\n")
  }
  
  # Run metaMDS, get points and stress
  require(vegan)
  mds <- metaMDS(comm, distance = distance)
  pts.mds <- as.data.frame(mds$points)
  pts.mds <- pts.mds[order(rownames(pts.mds)),]
  
  if (title != "")
    title <- paste0(title, " (stress ", round(mds$stress, 2), ")")
  
  if (! missing(attr.df)) {
    #rownames(pts.mds) <- tolower(rownames(pts.mds))
    #rownames(attr.df) <- tolower(rownames(attr.df))
    
    if (! is.null(colour.id)) {
      if (! colour.id %in% colnames(attr.df))
        stop("Invalid colour.id,", colour.id,  "not exsit in column names !\n")
    } 
    if (! is.null(link.id)) { 
      if (! link.id %in% colnames(attr.df) )
        stop("Invalid link.id,", link.id,  "not exsit in column names !\n")
    } 
    if (! is.null(shape.id)) { 
      if (! shape.id %in% colnames(attr.df) )
        stop("Invalid shape.id,", shape.id,  "not exsit in column names !\n")
    }
    
    pts.mds.merge <- merge(pts.mds, attr.df, by = "row.names")
    
    if (nrow(pts.mds.merge) != nrow(pts.mds) || nrow(pts.mds.merge) != nrow(attr.df)) 
      warning(paste("Some data points are missing after merge ! nrow(pts.mds.merge) =", 
                    nrow(pts.mds.merge), ", nrow(pts.mds) =", nrow(pts.mds), ", nrow(attr.df) =", nrow(attr.df) ))
  } else {
    pts.mds.merge <- pts.mds
    pts.mds.merge[,"Row.names"] <- rownames(pts.mds) 
  }
  
  if (add.text)
    text.id="Row.names"
  
  # Plot MDS ordination
  gt <- ComMA::gtScatterPlot(pts.mds.merge, x.id="MDS1", y.id="MDS2", colour.id=colour.id, 
                             shape.id=shape.id, link.id=link.id, text.id=text.id, 
                             palette=palette, xintercept=xintercept, yintercept=yintercept,
                             text.data=text.data, text.size=text.size, title=title, ...)
  
  return(gt)
}


#' PCA Plot
#' 
#' Use \code{\link{prcomp}} to create Principal Components Analysis (PCA) plot.
#' It is extended from \code{\link{gtScatterPlot}}. 
#' 
#' @param comm Community data for NMDS plot. 
#' It is dissimilarities either as a \code{\link{dist}} structure 
#' or as a community data defined by \code{\link{vegan}} which is 
#' a transposed matrix of community matrix in \code{\link{ComMA}}.
#' Use \code{\link{transpose.df}} to rotate the data frame.
#' Detail to \code{\link{prcomp}}.
#' @param attr.df A data frame to define how NMDS plot is coloured, shaped, or linked.
#' where its row names have to match row names of community data \code{comm}.
#' @keywords graph
#' @export
#' @examples
#' 
#' pca.plot <- gtNMDSPlot(comm, env, colour.id="FishSpecies", shape.id="FeedingPattern", add.text=T)
#' require(grid)
#' grid.draw(pca.plot)
gtPCAPlot <- function(comm, attr.df, x.i=1, y.i=2, colour.id=NULL, shape.id=NULL, link.id=NULL, 
                      add.text=TRUE, text.data=NULL, text.size=3, palette=NULL,
                      xintercept=NULL, yintercept=NULL, title="PCA", verbose=TRUE, ...) {
  if (! missing(attr.df)) {
    if (! all(rownames(as.matrix(comm)) %in% rownames(attr.df)) )
      stop("Invalid attr.df,", paste(rownames(as.matrix(comm)), collapse = ","), 
           "should match", paste(rownames(attr.df), collapse = ","), "!\n")
  }
  
  # Run prcomp, get points 
  pca <- prcomp(comm, scale. = TRUE)
  pts.pca <- as.data.frame(pca$rotation)
  pts.pca <- pts.pca[order(rownames(pts.pca)),]
  
  if (x.i>=y.i || x.i < 1 || y.i > ncol(pts.pca))
    stop("Invalid x.i", x.i,  "or y.i", y.i, "for PCA dimension index !\n")
  
  if (! missing(attr.df)) {
    #rownames(pts.pca) <- tolower(rownames(pts.pca))
    #rownames(attr.df) <- tolower(rownames(attr.df))
    
    if (! is.null(colour.id)) {
      if (! colour.id %in% colnames(attr.df))
        stop("Invalid colour.id,", colour.id,  "not exsit in column names !\n")
    } 
    if (! is.null(link.id)) { 
      if (! link.id %in% colnames(attr.df) )
        stop("Invalid link.id,", link.id,  "not exsit in column names !\n")
    } 
    if (! is.null(shape.id)) { 
      if (! shape.id %in% colnames(attr.df) )
        stop("Invalid shape.id,", shape.id,  "not exsit in column names !\n")
    }
    
    pts.pca.merge <- merge(pts.pca, attr.df, by = "row.names")
    
    if (nrow(pts.pca.merge) != nrow(pts.pca) || nrow(pts.pca.merge) != nrow(attr.df)) 
      warning(paste("Some data points are missing after merge !", 
                    nrow(pts.pca.merge), "!=", nrow(pts.pca), "!=", nrow(attr.df) ))
  } else {
    pts.pca.merge <- pts.pca
    pts.pca.merge[,"Row.names"] <- rownames(pts.pca) 
  }
  
  if (add.text)
    text.id="Row.names"
  
  x.id <- paste0("PC", x.i)
  y.id <- paste0("PC", y.i)
  # Plot MDS ordination
  gt <- ComMA::gtScatterPlot(pts.pca.merge, x.id="PC1", y.id="PC2", colour.id=colour.id, 
                             shape.id=shape.id, link.id=link.id, text.id=text.id, 
                             palette=palette, xintercept=xintercept, yintercept=yintercept,
                             text.data=text.data, text.size=text.size, title=title, ...)
  
  return(gt)
}

#' Percentage bar chart coloured by groups, which is extended from \code{\link{ggBarChart}}. 
#' 
#' @param df.to.melt A data frame required to \code{\link{melt}} before making a percent bar chart. 
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
#' bar.chart$gg.plot
#' pdfGgplot(bar.chart$gg.plot, fig.path="taxa-percentage-bar.pdf", width=bar.chart$pdf.width, height=8) 
ggPercentageBarChart <- function(df.to.melt, melt.id, title="Percentage Bar Chart", x.lab="", y.lab="", 
                                 x.text.angle=0, autoWidth=TRUE) {
  if (!is.element(tolower(melt.id), tolower(colnames(df.to.melt))))
    stop(paste0("Data frame column names do NOT have \"", melt.id, "\" for melt function !"))
  
  require(reshape2)
  df.melt <- melt(df.to.melt, id=c(melt.id))
  #df.melt[,"variable"] <- factor(df.melt[,"variable"], levels = sort(unique(df.melt[,"variable"])))
  
  # move unclassified group to the last of legend 
  legend.ord <- as.character(sort(unique(df.to.melt[,melt.id]), decreasing = TRUE))
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
  
  p <- ComMA::ggBarChart(df.melt, x.id="variable", y.id="value", fill.id=melt.id, bar.pos="fill", x.text.angle=x.text.angle,
                         y.trans="per", title=title, x.lab=x.lab, y.lab=y.lab, palette=pale, legend.col=legend.col)
  
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
