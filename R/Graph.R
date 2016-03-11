# Graph
# Author: Walter Xie
# Accessed on 10 Mar 2016

#' Create a heat map using ggplot. 
#' 
#' @param df A data frame to \code{\link{melt}} and then make a heat map. 
#' For example,
#' \tabular{rrrr}{
#'   plot \tab 16s \tab 18s \tab ITS\cr
#'   CM30c39 \tab 2 \tab 1 \tab 3\cr
#'   CM30c44 \tab 10 \tab 26 \tab 15\cr
#'   Plot01 \tab 6 \tab 5 \tab 6 
#' } 
#' @param id.melt A column name as a \code{\link{factor}}.
#' @param fig.path The full path of image file.
#' @param title Graph title
#' @param x.lab, y.lab The label of x-axis or y-axis, such as plot names
#' @param low, high Refer to \pkg{ggplot2} \code{\link{scale_fill_gradient}}. Default to low="white", high="steelblue".
#' @param width, height Refer to \code{\link{pdf}}. Default to width=6, height=6.
#' @keywords graph
#' @export
#' @examples 
#' ranks.by.group <- data.frame(plot=c("Plot03","Plot02","Plot01"), `16s`=c(3,2,1), `18s`=c(1,2,3), ITS=c(2,1,3), check.names = F)
#' ranks.by.group
#' heatmapGgplot(df=ranks.by.group, id.melt="plot", fig.path="plot-prior-example-heatmap.pdf")
heatmapGgplot <- function(df, id.melt, fig.path, title="Heatmap", x.lab="", y.lab="", 
                             low="white", high="steelblue", width=6, height=6) {
  if (!is.element(id.melt, tolower(colnames(df))))
    stop(paste0("Data frame column names do NOT have \"", id.melt, "\" for melt function !"))
  
  require(reshape2)
  
  breaks.rank <- round(seq(1, nrow(df), length.out = 5), digits = 0)
  ranks.melt <- melt(df, id=c(id.melt))
  ranks.melt[,id.melt] <- factor(ranks.melt[,id.melt], levels=unique(ranks.melt[,id.melt]))
  
  require(ggplot2)
  # variable is all group names, such as "16S" or "FUNGI"
  # value is ranks for each group
  p <- ggplot(ranks.melt, aes_string(x="variable", y=id.melt)) + geom_tile(aes(fill=value)) + 
    scale_fill_gradient(na.value="transparent", low=low, high=high, name="rank", breaks=breaks.rank) +
    xlab(x.lab) + ylab(y.lab) + ggtitle(title) +
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
  
  pdf(fig.path, width=width, height=height)	
  print(p)
  invisible(dev.off()) 
}

#' Data ellipse on a scatter plot. 
#' 
#' Superimposed data ellipse on a \pkg{ggplot2} scatter plot.
#' @source \url{http://stackoverflow.com/questions/2397097/how-can-a-data-ellipse-be-superimposed-on-a-ggplot2-scatterplot}.
#' 
#' @param df.clusters The 3-column data frame for plotting. 
#' The first 2 columns are coordinates, and 3rd is cluster names. 
#' @param addLabel Add row names to scatter plot. Default to FALSE
#' @param point.size The size of points from \code{\link{point.size}}. Default to 3.
#' @param palette The palette to colour clusters using \code{\link{scale_colour_brewer}}. 
#' Default to 'Set1' (max 8 colours). Refer to \url{http://www.datavis.ca/sasmac/brewerpal.html}.
#' @keywords graph
#' @export
#' @examples 
#' df.clusters <- random2Clusters()
#' df.clusters
#' scatterPlotEllipse(df.clusters, fig.path="clusters-scatter-plot.pdf")
#' scatterPlotEllipse(df.clusters, fig.path="clusters-scatter-plot.pdf", addLabel=T)
scatterPlotEllipse <- function(df.clusters, fig.path, title="Clusters", point.size=3, palette="Set1", 
                               addLabel=FALSE, label.size = 3, hjust=-0.1, vjust = -0.2, alpha = 0.5, 
                               width=8, height=6) {
  if (ncol(df.clusters) < 3)
    stop("Data frame should have 3 columns: first 2 columns are coordinates, 3rd is cluster names !")
  
  colnames(df.clusters)[1:3] <- c("PC1", "PC2", "cluster")
  # df.clusters$species <- paste(sapply(strsplit(rownames(df.clusters), "_"), "[[", 1), sapply(strsplit(rownames(df.clusters), "_"), "[[", 2), sep=".")
  df.clusters$cluster <- factor(df.clusters$cluster, levels = sort(unique(df.clusters$cluster)))
  
  require(ggplot2)
  
  p <- ggplot(df.clusters, aes(x=PC1, y=PC2, color=factor(cluster))) + 
    geom_point(size=point.size) + #aes(shape=factor(species)),
    scale_colour_brewer(name="cluster", palette=palette) + #scale_fill_manual(name = "cluster", values = myPalette)  +
    #    scale_shape_manual(name ="species", values=1:length(unique(gg$species))) +
    stat_ellipse(type = "t", linetype = 2) +
    geom_hline(yintercept=0,linetype=2) + 
    geom_vline(xintercept=0,linetype=2) +
    ggtitle(title) +
    theme_bw() + theme(panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), panel.background = element_blank()) 
  
  if (addLabel) {
    p <- p + geom_text(aes(color=factor(cluster), label=rownames(df.clusters)), size=label.size, hjust=hjust, vjust=vjust, alpha=alpha) 
  }
  
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  
  require(grid)
  
  pdf(fig.path, width=width, height=height)
  print(grid.draw(gt))
  invisible(dev.off()) 
}


