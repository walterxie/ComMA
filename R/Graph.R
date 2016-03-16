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
#' @param id.melt A column name to \code{\link{melt}} and used as a \code{\link{factor}}.
#' @param title Graph title
#' @param x.lab, y.lab The label of x-axis or y-axis, such as plot names.
#' @param low, high Refer to \pkg{ggplot2} \code{\link{scale_fill_gradient}}. Default to low="white", high="steelblue".
#' @return
#' A \code{\link{ggplot}} object of heatmap.
#' @keywords graph
#' @export
#' @examples 
#' ranks.by.group <- data.frame(plot=c("Plot03","Plot02","Plot01"), `16s`=c(3,2,1), `18s`=c(1,2,3), ITS=c(2,1,3), check.names = F)
#' ranks.by.group
#' gg.plot <- ggHeatmap(df=ranks.by.group, id.melt="plot")
#' pdfGgplot(gg.plot, fig.path="plot-prior-example-heatmap.pdf") 
ggHeatmap <- function(df, id.melt, title="Heatmap", x.lab="", y.lab="", low="white", high="steelblue") {
  if (!is.element(tolower(id.melt), tolower(colnames(df))))
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
  
  return(p) 
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
#' @return
#' A \code{\link{gtable}} object of scatter plot.
#' @keywords graph
#' @export
#' @examples 
#' df.clusters <- random2Clusters()
#' df.clusters
#' g.table <- ggScatterPlotEllipse(df.clusters, addLabel=T)
#' pdfGtable(g.table, fig.path="clusters-scatter-plot.pdf") 
ggScatterPlotEllipse <- function(df.clusters, title="Clusters", point.size=3, palette="Set1", 
                               addLabel=FALSE, label.size = 3, hjust=-0.1, vjust = -0.2, alpha = 0.5) {
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
  
  return(gt)
}

#' Percent bar chart coloured by groups. 
#' 
#' @param df A data frame to \code{\link{melt}} and then make a percent bar chart. 
#' For example,
#' \tabular{rrrrr}{
#'   Phyla \tab 16s \tab 18s \tab ITS \tab TaxaGroup\cr
#'   Actinobacteria \tab 958 \tab 1 \tab 3 \tab Bacteria\cr
#'   Crenarchaeota \tab 1 \tab 0 \tab 0 \tab Archaea\cr
#'   Ascomycota \tab 2 \tab 765 \tab 971 \tab Fungi 
#' } 
#' @param id.melt A column name to \code{\link{melt}} and used to assign the colours.
#' @param fig.path The full path of image file.
#' @param title Graph title
#' @param x.lab, y.lab The label of x-axis or y-axis, such as plot names.
#' @param low, high Refer to \pkg{ggplot2} \code{\link{scale_fill_gradient}}. 
#' Default to low="white", high="steelblue".
#' @param autoWidth If TRUE, then use number of bars and legend columns 
#' to estimate pdf width automatically. Default to TRUE.
#' @keywords graph
#' @export
#' @examples 
#' taxa.phyla <- readFile("./data/examples/taxonomy97phyla.txt")
#' bar.chart <- ggPercentBarChart(taxa.phyla, id.melt="TaxaGroup")
#' pdfGgplot(bar.chart$gg.plot, fig.path="taxa-percentage-bar.pdf", width=bar.chart$pdf.width, height=8) 
ggPercentBarChart <- function(df, id.melt, fig.path, title="Percent Bar Chart", x.lab="", y.lab="", autoWidth=TRUE) {
  if (!is.element(tolower(id.melt), tolower(colnames(df))))
    stop(paste0("Data frame column names do NOT have \"", id.melt, "\" for melt function !"))
  
  require(reshape2)
  df.melt <- melt(df, id=c(id.melt))
  #df.melt[,"variable"] <- factor(df.melt[,"variable"], levels = sort(unique(df.melt[,"variable"])))
  
  # move unclassified group to the last of legend 
  legend.ord <- as.character(unique(df[,id.melt]))
  id.match <- grep("unclassified", legend.ord, ignore.case = TRUE)
  if (length(id.match) > 0)
    legend.ord <- legend.ord[c(setdiff(1:length(legend.ord), id.match),id.match)]
  df.melt[,id.melt] <- factor(df.melt[,id.melt], levels = rev(legend.ord))
  
  pale <- ComMA::getMyPalette(length(legend.ord))
  if (length(legend.ord) > length(pale)) {
    require(colorspace)
    pale <- rainbow_hcl(length(legend.ord))
  }
  # number of columns for legend
  legend.col = ceiling(length(legend.ord) / 25)
  
  require(ggplot2)
  require(scales)
  p <- ggplot(df.melt, aes_string(x = "variable", y = "value", fill = id.melt)) + 
    geom_bar(position = "fill",stat = "identity") + 
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(values=pale) +
    guides(fill=guide_legend(ncol=legend.col)) +
    theme_bw() + xlab(x.lab) + ylab(y.lab) + ggtitle(title) +
    theme(axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank()) 
  
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


#' Grouping bar chart Y across X. 
#' 
#' Create a grouping bar chart given community matrix to display 
#' the number of OTUs (y-axis) across the number of samples (x-axis). 
#' The 'red' bar is the number of OTUs appeared only in that number of samlpes, 
#' and the 'green' bar is the number of reads assigned to those OTUs.
#' @param df.aggre A data frame of row counts and sums from a community matrix. See \code{\link{rowCountAndSum}}. 
#' @param print.xtable TRUE/FALSE to print \code{\link{xtable}} in console. 
#' Allow NULL, if NULL, then do not print. Default to NULL. 
#' @param title Graph title
#' @param x.lab, y.lab The label of x-axis or y-axis, such as plot names.
#' @param legend.title The title of legend. Refer to \pkg{ggplot2} \code{\link{scale_fill_discrete}}. 
#' Default to a empty string.
#' @param legend.labels The labels of legend, which are fixed to 2 groups. Default to c("OTUs", "reads").
#' @param x.lab.interval The x labels interval. Default to 1 to dispay lables for every values.
#' @keywords graph
#' @export
#' @examples 
#' communityMatrix <- getCommunityMatrix("16S", isPlot=TRUE, minAbund=1)
#' cm.aggre <- cmYAcrossX(communityMatrix)
#' gg.plot <- ggBarYAcrossX(cm.aggre)
#' pdfGgplot(gg.plot, fig.path="bar-otus-across-samples.pdf", width=8, height=8)  
ggBarYAcrossX <- function(df.aggre, id.melt="samples", print.xtable=NULL, title="The number of OTUs across the number of samples", 
                       x.lab="Number of samples crossed", y.lab="Number of OTUs", 
                       legend.title="", legend.labels=c("OTUs", "reads"), x.lab.interval=1) {
  if (!is.element(tolower(id.melt), tolower(colnames(df.aggre))))
    stop(paste0("Data frame column names do NOT have \"", id.melt, "\" for melt function !"))
  
  require(reshape2)
  df.melt <- melt(df.aggre, id=id.melt)
  
  x.breaks <- seq(min(df.aggre[,id.melt]), max(df.aggre[,id.melt]), x.lab.interval)
  y.breaks <- ComMA::get_breaks_positive_values(max(df.aggre, start=c(0)))
  require(ggplot2)
  # conside x as discrete values
  p <- ggplot(df.melt, aes_string(x=id.melt, y="value")) + 
    geom_bar(aes(fill=variable), position = "dodge", stat="identity") +
    scale_y_continuous(trans="log", expand = c(0,0), labels = ComMA::scientific_10, breaks = y.breaks) + 
    scale_x_discrete(breaks=x.breaks) +
    xlab(x.lab) + ylab(y.lab) + ggtitle(title) +
    scale_fill_discrete(legend.title, labels=legend.labels) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black", size = 0.3),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  return(p)
}


