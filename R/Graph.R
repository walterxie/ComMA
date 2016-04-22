# Graph
# Author: Walter Xie
# Accessed on 12 Apr 2016

#' @name ggPlot
#' @title One-line command to get \pkg{ggplot} 
#'
#' @description Simplify \pkg{ggplot} codes into functions that can get a chart from one-line command.
#'
#' @return 
#' If the function returns a \code{\link{ggplot}} object, then its name starts with "gg". 
#' It needs to use \code{\link{pdfGgplot}} to create pdf. 
#' If the function returns a \code{\link{gtable}} object, then its name starts with "gt".
#' It needs to use \code{\link{pdfGtplot}} to create pdf. 
#' 

#' @details 
#' \code{ggAddLine} adds a line to the given \code{\link{ggplot}} object.
#' 
#' @param xintercept,yintercept,intercept,slope,smooth.method Refer to \pkg{ggplot2} 
#' \code{\link{geom_vline}}, \code{\link{geom_hline}}, \code{\link{geom_abline}}, 
#' \code{\link{geom_smooth}}. They cannot be used at the same time.
#' @param linetype \code{\link{linetype}} 0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash.
#' @keywords graph
#' @export
#' @examples 
#' p <- ggAddLine(p, linetype = 2, yintercept = 1)
#' p <- ggAddLine(p, smooth.method = "lm")
#' 
#' @rdname ggPlot
ggAddLine <- function(gg.plot, linetype=1, xintercept, yintercept, intercept, slope, smooth.method) {
  if (!missing(xintercept)) {
    p <- gg.plot + geom_vline(linetype=linetype, xintercept = xintercept)
  } else if (!missing(yintercept)) {
    p <- gg.plot + geom_hline(linetype=linetype, yintercept = yintercept)
  } else if (!missing(intercept) && !missing(slope)){
    p <- gg.plot + geom_abline(linetype=linetype, intercept = intercept, slope = slope)
  } else if (!missing(smooth.method)){  
    p <- p + geom_smooth(linetype=linetype, method = smooth.method, se = FALSE)
  } else {
    stop("Invalid input !")
  }
}

#' @details 
#' \code{ggAddNumbers} adds numbers as text in a \code{\link{ggplot}} object, such as mean of box plot. 
#' Refer to \code{\link{stat_summary}}.
#' 
#' @param fun.y.lab A function to calculate numbers displayed in the figure.  
#' Default to function \code{\link{mean}}. Ues \code{\link{length}} to show number of observations.
#' @param fun.y.pos A function to calculate the initial poistion of text on y-value. 
#' Default to \code{\link{median}}.
#' @param y.adj The propotion of the initial poistion of text on y-value. 
#' > 1 will raises the text, and < 1 will sinks the text. Default to 0.98.
#' @param digits Integer indicating the number of decimal places for \code{\link{round}}.
#' @param dodge.width Dodging width, when different to the width of the individual elements. 
#' Default to 0.8. Refer to \code{\link{position_dodge}}.
#' @param text.size The text size of labels. Default to 3.
#' @param text.colour The text colour. Default to black.
#' @keywords graph
#' @export
#' @examples 
#' p <- ggAddNumbers(p, fun.y.lab=mean)
#' p <- ggAddNumbers(p, fun.y.lab=length, y.adj=1.02)
#' 
#' @rdname ggPlot
ggAddNumbers <- function(gg.plot, fun.y.lab=mean, fun.y.pos=median, y.adj=0.98, digits=2, 
                         dodge.width=0.8, text.size=3, text.colour="black") {
  p <- gg.plot + stat_summary(fun.data = function(y) {return( c(y = fun.y.pos(y)*y.adj, label = round(fun.y.lab(y),digits)) )}, 
                              geom = "text", position = position_dodge(width=dodge.width), colour = text.colour, size = text.size)
}

#' @details 
#' \code{ggHeatmap} creates a heat map using ggplot. 
#' 
#' @param df A data frame to \code{\link{melt}} and then make a heat map. 
#' For example,
#' \tabular{rrrr}{
#'   plot \tab 16s \tab 18s \tab ITS\cr
#'   CM30c39 \tab 2 \tab 1 \tab 3\cr
#'   CM30c44 \tab 10 \tab 26 \tab 15\cr
#'   Plot01 \tab 6 \tab 5 \tab 6 
#' } 
#' @param melt.id A column name to \code{\link{melt}} and used as a \code{\link{factor}}.
#' @param title Graph title
#' @param x.lab,y.lab The label of x-axis or y-axis, such as plot names.
#' @param low,high Refer to \pkg{ggplot2} \code{\link{scale_fill_gradient}}. Default to low="white", high="steelblue".
#' @return
#' A \code{\link{ggplot}} object of heatmap.
#' @keywords graph
#' @export
#' @examples 
#' ranks.by.group <- data.frame(plot=c("Plot03","Plot02","Plot01"), `16s`=c(3,2,1), `18s`=c(1,2,3), ITS=c(2,1,3), check.names = F)
#' ranks.by.group
#' gg.plot <- ggHeatmap(df=ranks.by.group, melt.id="plot")
#' pdfGgplot(gg.plot, fig.path="plot-prior-example-heatmap.pdf") 
#' 
#' @rdname ggPlot
ggHeatmap <- function(df, melt.id, title="Heatmap", x.lab="", y.lab="", low="white", high="steelblue") {
  if (!is.element(tolower(melt.id), tolower(colnames(df))))
    stop("Data frame column names do NOT have \"", melt.id, "\" for melt function !")
  
  require(reshape2)
  breaks.rank <- round(seq(1, nrow(df), length.out = 5), digits = 0)
  ranks.melt <- melt(df, id=c(melt.id))
  ranks.melt[,melt.id] <- factor(ranks.melt[,melt.id], levels=unique(ranks.melt[,melt.id]))
  
  require(ggplot2)
  # variable is all group names, such as "16S" or "FUNGI"
  # value is ranks for each group
  p <- ggplot(ranks.melt, aes_string(x="variable", y=melt.id)) + geom_tile(aes(fill=value)) + 
    scale_fill_gradient(na.value="transparent", low=low, high=high, name="rank", breaks=breaks.rank) +
    xlab(x.lab) + ylab(y.lab) + ggtitle(title) +
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
  
  return(p) 
}

#' @details 
#' \code{ggBarChart} is an one-line function to plot many types of bar chart, 
#' such as normal bars, log-scaled bars, percentage bars, and also grouping.
#' 
#' @param df.melt A data frame already \code{\link{melt}}. 
#' @param x.id,y.id,fill.id,group.id The string of column names in \code{df.melt},
#' which use for \code{x, y, fill, group} in \code{\link{aes}} in \code{\link{ggplot}}.
#' @param x.facet.id, y.facet.id The string of column names in \code{df.melt},
#' which creates facets (a formula) in \code{\link{facet_grid}} in \code{\link{ggplot}}.
#' @param bar.pos Position adjustment for \code{\link{geom_bar}}, either as a string, 
#' or the result of a call to a position adjustment function. Default to "dodge". 
#' Use \code{fill} to generate group percentage bars.
#' @param bar.stat Determine what is mapped to bar height. Refer to \code{\link{geom_bar}}. 
#' Default to "identity", which define the heights of the bars to represent values in the data.
#' @param x.trans,y.trans The string defines the data scale used in either x-axis or y-axis, 
#' which can be "identity" standing for normal, or "per" standing for percentage, 
#' moreover either the name of a transformation object for \code{\link{scale_x_continuous}}
#' or \code{\link{scale_y_continuous}} (e.g. \code{trans="log"}), or the object itself. 
#' Built-in transformations include "asn", "atanh", "boxcox", "exp", "identity", 
#' "log", "log10", "log1p", "log2", "logit", "probability", "probit", "reciprocal", 
#' "reverse" and "sqrt". Default to "identity". 
#' @param x.lim.cart,y.lim.cart Setting limits on the coordinate system will zoom the plot, 
#' and will not change the underlying data like setting limits on a scale will. 
#' Refer to \code{\link{coord_cartesian}}. Set lower bound only to y-axis using y.lim.cart=c(1000,NA). 
#' Default NULL. 
#' @param title Graph title, set title="" to remove it from the plot.
#' @param x.lab,y.lab The label of x-axis or y-axis, such as plot names. 
#' Set x.lab="" to remove x-axis label from the plot.
#' @param x.lab.interval The interval to display x values in axis. 
#' Assume x values are discrete for each bar. Default to 0 to do nothing.
#' @param legend.title The title of legend. Set legend.title="" to remove legend.
#' @param palette The colour palette for bar, box, scatter plot, etc. 
#' If length == 1, then use \code{\link{scale_colour_brewer}} 
#' (\url{http://www.datavis.ca/sasmac/brewerpal.html}), such as "Set1" (max 8 colours).
#' If 1 != length <= 3, then use \code{\link{scale_colour_gradientn}}, 
#' such as c("blue", "orange").
#' Otherwise use \code{\link{scale_fill_manual}} for a vector of customized colours.
#' Default NULL to use \code{\link{ggplot}} default colours.  
#' @param legend.col,legend.nrow Customize the number of columns or rows for legend in bar chart. 
#' They cannot be used at the same time. Default not to use them, legend.col=1, legend.nrow=0. 
#' @keywords graph
#' @export
#' @examples
#' # log-scale y
#' bar.chart <- ggBarChart(df.melt, x.id="test", y.id="seconds", fill.id="version", y.trans="log")
#' # percentage bars without grouping in one bar each
#' bar.chart <- ggBarChart(df.melt, x.id="test", y.id="percentage", fill.id="model", y.trans="per")
#' # percentage bars one group in one bar
#' bar.chart <- ggBarChart(df.melt, x.id="test", y.id="percentage", fill.id="model", bar.pos="fill", y.trans="per")
#' 
#' # the number of OTUs (y-axis) across the number of samples (x-axis)
#' communityMatrix <- getCommunityMatrix("16S", isPlot=TRUE, minAbund=1)
#' cm.aggre <- cmYAcrossX(communityMatrix)
#' require(reshape2)
#' df.melt <- melt(cm.aggre, id="samples")
#' bar.chart <- ggBarChart(df.melt, x.id="samples", y.id="value", fill.id="variable", y.trans="log", 
#'                         y.lab="", legend.title="", x.interval=1)
#' 
#' @rdname ggPlot
ggBarChart <- function(df.melt, x.id, y.id, fill.id=NULL, bar.pos="dodge", bar.stat="identity", 
                       x.facet.id=NULL, y.facet.id=NULL, y.trans="identity", auto.scale.y=FALSE, 
                       x.interval=0, x.lim.cart=NULL, y.lim.cart=NULL, palette=NULL, 
                       legend.title=NULL, legend.col=1, legend.nrow=0, x.text.angle=0, 
                       title="Bar Chart", title.size = 10, x.lab="x.id", y.lab="y.id") {
  p <- ggInit(df.melt=df.melt, x.id=x.id, y.id=y.id, fill.id=fill.id)
  p <- p + geom_bar(position=bar.pos, stat=bar.stat) 
  
  col.names <- colnames(df.melt)
  p <- ggOptFacetGrid(p, col.names, x.facet.id=x.facet.id, y.facet.id=y.facet.id)
  
  if (auto.scale.y) {
    y.max <- max(df.melt[,y.id])
    p <- ggOptScaleAxis(p, axis="y", scale="continuous", trans=y.trans, auto.scale.max=y.max)
  } else {
    p <- ggOptScaleAxis(p, axis="y", scale="continuous", trans=y.trans)
  }
  
  if (x.interval > 0) {
    x.breaks <- seq(min(df.melt[,x.id]), max(df.melt[,x.id]), x.interval)
    p <- ggOptScaleAxis(p, axis="x", scale="discrete", breaks=x.breaks)
  }
  
  p <- ggOptCoordCartesian(p, df.melt, x.id, y.id, x.lim.cart=x.lim.cart, y.lim.cart=y.lim.cart)
  
  p <- ggOptPalette(p, palette=palette)
  
  p <- ggOptLegend(p, legend.title=legend.title, legend.col=legend.col, legend.nrow=legend.nrow)
  
  p <- ggThemePanelBorder(p, x.id, y.id, title=title, title.size=title.size, x.lab=x.lab, y.lab=y.lab)
  
  p <- ggThemeRotateXText(p, x.text.angle=x.text.angle)
  
  return(p)
}


#' @details 
#' \code{ggBoxWhiskersPlot} creates box Whiskers plot. 
#' 
#' @param outlier.colour The colour of outliers in box whiskers plot 
#' used for \code{outlier.colour} in \code{\link{geom_boxplot}} in \code{\link{ggplot}}. 
#' Default to alpha("black", 0.3).
#' @param dodge.width Dodging width, when different to the width of the individual elements. 
#' This is useful when you want to align narrow geoms with wider geoms. 
#' Refer to \code{\link{position_dodge}}.
#' @keywords graph
#' @export
#' @examples
#' box.plot <- ggBoxWhiskersPlot(df.melt, x.id="test", y.id="performance")
#' 
#' @rdname ggPlot
ggBoxWhiskersPlot <- function(df.melt, x.id, y.id, fill.id=NULL, 
                              outlier.colour=alpha("black", 0.3), dodge.width=0.8,
                              x.facet.id=NULL, y.facet.id=NULL, 
                              y.trans="identity", auto.scale.y=FALSE, 
                              x.lim.cart=NULL, y.lim.cart=NULL, palette=NULL, 
                              legend.title=NULL, legend.col=1, legend.nrow=0, x.text.angle=0, 
                              title="Box Whiskers Plot", title.size = 10, x.lab="x.id", y.lab="y.id") {
  p <- ggInit(df.melt=df.melt, x.id=x.id, y.id=y.id, fill.id=fill.id)
  if (! is.null(fill.id)) 
    p <- p + geom_boxplot(outlier.colour=outlier.colour, position=position_dodge(width=dodge.width))
  else 
    p <- p + geom_boxplot(outlier.colour=outlier.colour)
  
  p <- p + scale_shape(solid = FALSE) #+ geom_jitter(alpha = 0.5) 
  
  col.names <- colnames(df.melt)
  p <- ggOptFacetGrid(p, col.names, x.facet.id=x.facet.id, y.facet.id=y.facet.id)
  
  if (auto.scale.y) {
    y.max <- max(df.melt[,y.id])
    p <- ggOptScaleAxis(p, axis="y", scale="continuous", trans=y.trans, auto.scale.max=y.max)
  } else {
    p <- ggOptScaleAxis(p, axis="y", scale="continuous", trans=y.trans)
  }
  
  p <- ggOptCoordCartesian(p, df.melt, x.id, y.id, x.lim.cart=x.lim.cart, y.lim.cart=y.lim.cart)
  
  p <- ggOptPalette(p, palette=palette)
  
  p <- ggOptLegend(p, legend.title=legend.title, legend.col=legend.col, legend.nrow=legend.nrow)
  
  p <- ggThemePanelBorder(p, x.id, y.id, title=title, title.size=title.size, x.lab=x.lab, y.lab=y.lab)
  
  p <- ggThemeRotateXText(p, x.text.angle=x.text.angle)
  
  return(p)
}


#' @details 
#' \code{gtScatterPlot} uses one-line function to plot many types of scatter chart.
#' 
#' @param point.size The size of points from \code{\link{point.size}}. Default to 3.
#' @param text.id Label the data points according \code{text.id} column, 
#' such as "Row.names" column after \code{\link{merge}}.
#' @param text.size,text.hjust,text.vjust,text.alpha 
#' The parameters to adjust text in \code{\link{geom_text}}.
#' @param colour.id,shape.id,linke.id The column name in \code{df.melt} to 
#' define how the data points are coloured, shaped, or linked according their values.
#' @param ellipsed.id The column name in \code{df.melt} to define 
#' how to draw ellipse over data points, which is normally same as 
#' \code{colour.id} to show clusters.
#' @param shapes Manually define the shapes of points. Refer to \code{\link{scale_shape_manual}}, 
#' and \url{http://www.cookbook-r.com/Graphs/Shapes_and_line_types/}.
#' @param xintercept,yintercept,linetype Add horizontal or vertical line. 
#' Refer to \code{\link{geom_hline}} or \code{\link{geom_vline}}.
#' @keywords graph
#' @export
#' @examples 
#' df.clusters <- random2Clusters()
#' df.clusters$labels <- rownames(df.clusters)
#' df.clusters
#' g.table <- gtScatterPlot(df.clusters, x.id="x", y.id="y", colour.id="group", shape.id="group",   
#'                          xintercept=0, yintercept=0, title="Clusters", palette="Set1")
#' require(grid)
#' grid.draw(g.table)
#' # selective labeling for points x > 3 and y > 6
#' g.table <- gtScatterPlot(df.clusters, x.id="x", y.id="y", colour.id="group", ellipsed.id="group",
#'                          text.id="labels", text.data=subset(df.clusters, x > 3 & y > 6), 
#'                          xintercept=0, yintercept=0, title="Clusters", palette="Set1")
#' grid.draw(g.table)
#' pdfGtable(g.table, fig.path="clusters-scatter-plot.pdf") 
#' 
#' @rdname ggPlot
gtScatterPlot <- function(df.melt, x.id, y.id, colour.id=NULL, shape.id=NULL, 
                          shapes=NULL, point.size=3, x.facet.id=NULL, y.facet.id=NULL, 
                          linke.id=NULL, ellipsed.id=NULL, text.id=NULL, 
                          text.data = NULL, text.size = 3, text.hjust=-0.1, 
                          text.vjust = -0.2, text.alpha = 0.5,
                          xintercept=NULL, yintercept=NULL, line.type=2,
                          x.lim.cart=NULL, y.lim.cart=NULL, palette=NULL,
                          legend.title=NULL, legend.col=1, legend.nrow=0, 
                          title="Scatter Plot", title.size = 10, x.lab="x.id", y.lab="y.id") {
  p <- ggInit(df.melt=df.melt, x.id=x.id, y.id=y.id, colour.id=colour.id)
  col.names <- colnames(df.melt)
  
  if (! is.null(shape.id) && is.null(shapes)) {
    # The shape palette can deal with a maximum of 6 discrete values
    n_shape <- length(unique(df.melt[,shape.id]))
    shapes <- seq(1, (1 + n_shape-1))
  }
  p <- ggOptPointAndShape(p, col.names, shape.id=shape.id, shapes=shapes, point.size=point.size)
  
  p <- ggOptFacetGrid(p, col.names, x.facet.id=x.facet.id, y.facet.id=y.facet.id)
  
  p <- ggOptEllipse(p, col.names, ellipsed.id=ellipsed.id)
  
  if (! is.null(linke.id)) {
    # Convex hull http://stackoverflow.com/questions/16428962/convex-hull-ggplot-using-data-tables-in-r
    df.melt <- data.table(df.melt, key = linke.id)
    chull.cmd <- parse(text = paste0('df.melt[, .SD[chull(', x.id, ',', y.id, ')], by = "', linke.id, '"')) 
    # hulls <- pts.mds.dt[, .SD[chull(MDS1, MDS2)], by = linke.id]
    hulls <- eval(chull.cmd)
    
    p <- p + geom_polygon(data = hulls, aes_string(mapping=linke.id), fill = NA, alpha = 0.5)
  }
  
  p <- ggOptText(p, col.names, text.id=text.id, text.data=text.data, colour.id=colour.id, text.size=text.size, 
                 text.hjust=text.hjust, text.vjust=text.vjust, text.alpha=text.alpha)
  
  if (! is.null(xintercept))
    p <- p + geom_vline(xintercept=xintercept,linetype=line.type)
  if (! is.null(yintercept))
    p <- p + geom_hline(yintercept=yintercept,linetype=line.type) 
  
  p <- ggOptCoordCartesian(p, df.melt, x.id, y.id, x.lim.cart=x.lim.cart, y.lim.cart=y.lim.cart)
  
  p <- ggOptPalette(p, palette=palette)
  
  p <- ggOptLegend(p, legend.title=legend.title, legend.col=legend.col, legend.nrow=legend.nrow)
  
  p <- ggThemePanelBorder(p, x.id, y.id, title=title, title.size=title.size, x.lab=x.lab, y.lab=y.lab)
  
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  
  return(gt)
}

#' @details 
#' \code{gtLine} uses one-line function to plot a line or group of lines.
#' 
#' @param line.size,line.alpha The feature of lines for \code{\link{geom_line}}.
#' @keywords graph
#' @export
#' @examples 
#' 
#' 
#' @rdname ggPlot
gtLine <- function(df.melt, x.id, y.id, group.id=NULL, colour.id=NULL,  
                   line.size=0.5, line.type = 2, line.alpha=0.75, 
                   shape.id=NULL, shapes=NULL, point.size=3, 
                   x.facet.id=NULL, y.facet.id=NULL, y.trans="identity", auto.scale.y=FALSE,
                   text.id=NULL, text.data = NULL, text.size = 3, 
                   text.hjust=-0.1, text.vjust = -0.2, text.alpha = 0.5, 
                   x.lim.cart=NULL, y.lim.cart=NULL, palette=NULL, 
                   legend.title=NULL, legend.col=1, legend.nrow=0, x.text.angle=0, 
                   title="", title.size = 10, x.lab="x.id", y.lab="y.id") {
  p <- ggInit(df.melt=df.melt, x.id=x.id, y.id=y.id, group.id=group.id, colour.id=colour.id)
  col.names <- colnames(df.melt)
  
  p <- p + geom_line(size=line.size, linetype=line.type, alpha=line.alpha) 
  
  if (! is.null(shape.id) && is.null(shapes)) {
    # The shape palette can deal with a maximum of 6 discrete values
    n_shape <- length(unique(df.melt[,shape.id]))
    shapes <- seq(1, (1 + n_shape-1))
  }
  p <- ggOptPointAndShape(p, col.names, shape.id=shape.id, shapes=shapes, point.size=point.size)
  
  if (auto.scale.y) {
    y.max <- max(df.melt[,y.id])
    p <- ggOptScaleAxis(p, axis="y", scale="continuous", trans=y.trans, auto.scale.max=y.max)
  } else {
    p <- ggOptScaleAxis(p, axis="y", scale="continuous", trans=y.trans)
  }
  
  p <- ggOptText(p, col.names, text.id=text.id, text.data=text.data, colour.id=colour.id, text.size=text.size, 
                 text.hjust=text.hjust, text.vjust=text.vjust, text.alpha=text.alpha)
  
  p <- ggOptCoordCartesian(p, df.melt, x.id, y.id, x.lim.cart=x.lim.cart, y.lim.cart=y.lim.cart)
  
  p <- ggOptPalette(p, palette=palette)
  
  p <- ggOptLegend(p, legend.title=legend.title, legend.col=legend.col, legend.nrow=legend.nrow)
  
  p <- ggThemeAxis(p, x.id, y.id, title=title, title.size=title.size, x.lab=x.lab, y.lab=y.lab)
  
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  
  return(gt)
}

