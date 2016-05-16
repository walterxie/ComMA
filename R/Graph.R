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
#' It needs to use \code{\link{pdf.ggplot}} to create pdf. 
#' If the function returns a \code{\link{gtable}} object, then its name starts with "gt".
#' It needs to use \code{\link{pdf.gtplot}} to create pdf. 
#' 
#' @param df A data frame used for plot. 
#' @param df.to.melt A data frame required to \code{\link{melt}} 
#' before making a \pkg{ggplot} object, such as input of \code{ggHeatmap}. 
#' At least one column should be \code{melt.id}. If using row.names, 
#' then it should be inserted into data frame before this function. 
#' For example,
#' \tabular{rrrr}{
#'   plot \tab 16s \tab 18s \tab ITS\cr
#'   CM30c39 \tab 2 \tab 1 \tab 3\cr
#'   CM30c44 \tab 10 \tab 26 \tab 15\cr
#'   Plot01 \tab 6 \tab 5 \tab 6 
#' } 
#' @param melt.id A column name to \code{\link{melt}} 
#' and used as a \code{\link{factor}}, such as "plot" column.
#' @param x.id,y.id,fill.id,group.id The string of column names in \code{df},
#' which use for \code{x, y, fill, group} in \code{\link{aes}} in \code{\link{ggplot}}.
#' @param x.facet.id, y.facet.id The string of column names in \code{df},
#' which creates facets (a formula) in \code{\link{facet_grid}} in \code{\link{ggplot}}.
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
#' @param coord.flip If TRUE, then flip cartesian coordinates so that horizontal 
#' becomes vertical, and vertical becomes horizontal. Default to FALSE. Refer to 
#' \code{\link{coord_flip}}.
#' @param x.lab.interval The interval to display x values in axis. 
#' Assume x values are discrete for each bar. Default to 0 to do nothing.
#' @param palette The colour palette for bar, box, scatter plot, etc. 
#' If length == 1, then use \code{\link{scale_colour_brewer}} 
#' (\url{http://www.datavis.ca/sasmac/brewerpal.html}), such as "Set1" (max 8 colours).
#' If 1 != length <= 3, then use \code{\link{scale_colour_gradientn}}, 
#' such as c("blue", "orange").
#' Otherwise use \code{\link{scale_fill_manual}} for a vector of customized colours.
#' Default NULL to use \code{\link{ggplot}} default colours.  
#' @param text.id Label the data points according \code{text.id} column, 
#' such as "Row.names" column after \code{\link{merge}}.
#' @param text.size,text.hjust,text.vjust,text.alpha 
#' The parameters to adjust text in \code{\link{geom_text}}.
#' @param legend.title The title of legend. Set legend.title="" to remove legend.
#' @param legend.col,legend.row Customize the number of columns or rows for legend in bar chart. 
#' They cannot be used at the same time. Default not to use them, legend.col=1, legend.row=0. 
#' @param no.panel.border Add panel border or not. Default to FALSE.




#' @details 
#' \code{ggHeatmap} creates a heat map using ggplot. 
#' 
#' @param low,high Refer to \pkg{ggplot2} \code{\link{scale_fill_gradient}}. 
#' Default to low="white", high="steelblue".
#' @param log.scale.colour If TRUE, then use log scale to the colour of heat map.
#' Default to FALSE.
#' @param x.text,y.text If FALSE, then hide x or y axis labels. Default to TRUE.
#' @keywords graph
#' @export
#' @examples 
#' ranks.by.group <- data.frame(plot=c("Plot03","Plot02","Plot01"), `16s`=c(3,2,1), 
#'                              `18s`=c(1,2,3), ITS=c(2,1,3), check.names = F)
#' ranks.by.group
#' gg.plot <- ggHeatmap(ranks.by.group, melt.id="plot")
#' pdf.ggplot(gg.plot, fig.path="plot-prior-example-heatmap.pdf") 
#' 
#' @rdname ggPlot
ggHeatmap <- function(df.to.melt, melt.id, low="white", high="steelblue", 
                      title="Heatmap", title.size = 10, x.lab="", y.lab="", 
                      log.scale.colour=FALSE, legend.title="Counts",
                      x.lim.cart=NULL, y.lim.cart=NULL, coord.flip=FALSE,
                      x.text=TRUE, y.text=TRUE, x.text.angle=45, 
                      no.panel.border=FALSE, verbose=TRUE) {
  if (!is.element(tolower(melt.id), tolower(colnames(df.to.melt))))
    stop("Data frame column names do NOT have \"", melt.id, "\" for melt function !")
  
  suppressMessages(require(reshape2))
  df.melt <- melt(df.to.melt, id=c(melt.id))
  df.melt[,melt.id] <- factor(df.melt[,melt.id], levels=unique(df.melt[,melt.id]))
  
  suppressMessages(require(ggplot2))
  # variable is all group names, such as "16S" or "FUNGI"
  # value is ranks for each group
  p <- ggplot(df.melt, aes_string(x="variable", y=melt.id)) + geom_tile(aes(fill=value)) 
  
  if (log.scale.colour) {
    if (min(df.melt$value) == 0)
      min.log <- 0
    else 
      min.log <- log(min(df.melt$value))
    if (max(df.melt$value) == 0)
      max.log <- 0
    else 
      max.log <- log(max(df.melt$value))
    
    breaks.rank <- round(exp(seq(min.log, max.log, length.out = 5)), digits = 0)
    p <- p + scale_fill_gradient(trans='log', na.value="transparent", low=low, high=high, 
                                 name=legend.title, breaks=breaks.rank) 
  } else {
    breaks.rank <- round(seq(min(df.melt$value), max(df.melt$value), length.out = 5), digits = 0)
    p <- p + scale_fill_gradient(na.value="transparent", low=low, high=high, 
                                 name=legend.title, breaks=breaks.rank) 
  }
  
  p <- ggOptCoordCartesian(p, df.melt, x.id, y.id, x.lim.cart=x.lim.cart, y.lim.cart=y.lim.cart, 
                           coord.flip=coord.flip, verbose=verbose)
  
  p <- ggLabTitle(p, "", "", title=title, x.lab=x.lab, y.lab=y.lab)
  if (no.panel.border)
    p <- ggThemeAxis(p, title.size=title.size)
  else 
    p <- ggThemePanelBorder(p, title.size=title.size)
  
  p <- ggThemeOthers(p, x.text.angle=x.text.angle, x.text=x.text, y.text=y.text)
  
  return(p) 
}

#' @details 
#' \code{ggBarChart} is an one-line function to plot many types of bar chart, 
#' such as normal bars, log-scaled bars, percentage bars, and also grouping.
#' 
#' @param bar.pos Position adjustment for \code{\link{geom_bar}}, either as a string, 
#' or the result of a call to a position adjustment function. Default to "dodge". 
#' Use \code{fill} to generate group percentage bars.
#' @param bar.stat Determine what is mapped to bar height. Refer to \code{\link{geom_bar}}. 
#' Default to "identity", which define the heights of the bars to represent values in the data.
#' @keywords graph
#' @export
#' @examples
#' # log-scale y
#' bar.chart <- ggBarChart(df, x.id="test", y.id="seconds", fill.id="version", y.trans="log")
#' # percentage bars without grouping in one bar each
#' bar.chart <- ggBarChart(df, x.id="test", y.id="percentage", fill.id="model", y.trans="per")
#' # percentage bars one group in one bar
#' bar.chart <- ggBarChart(df, x.id="test", y.id="percentage", fill.id="model", bar.pos="fill", y.trans="per")
#' 
#' @rdname ggPlot
ggBarChart <- function(df, x.id, y.id, fill.id=NULL, bar.pos="dodge", bar.stat="identity", 
                       x.facet.id=NULL, y.facet.id=NULL, x.lim.cart=NULL, y.lim.cart=NULL,   
                       y.trans="identity", auto.scale.y=FALSE, x.scale="discrete", 
                       x.interval=0, x.text.angle=0, palette=NULL, coord.flip=FALSE,
                       legend.title=NULL, legend.col=1, legend.row=0, 
                       title="Bar Chart", title.size=10, x.lab="x.id", y.lab="y.id", 
                       legend.position="right", legend.direction="vertical",
                       no.panel.border=FALSE, verbose=TRUE) {
  p <- ggInit(df=df, x.id=x.id, y.id=y.id, fill.id=fill.id, verbose=verbose)
  p <- p + geom_bar(position=bar.pos, stat=bar.stat) 
  
  col.names <- colnames(df)
  p <- ggOptFacetGrid(p, col.names, x.facet.id=x.facet.id, y.facet.id=y.facet.id)
  
  if (auto.scale.y) {
    y.max <- max(df[,y.id])
    p <- ggOptScaleAxis(p, axis="y", scale="continuous", trans=y.trans, 
                        auto.scale.max=y.max, verbose=verbose)
  } else {
    p <- ggOptScaleAxis(p, axis="y", scale="continuous", trans=y.trans, 
                        verbose=verbose)
  }
  
  if (x.interval > 0) {
    #x.breaks <- seq(min(df[,x.id]), max(df[,x.id]), x.interval)
    x.breaks <- window(unique(df[,x.id]), deltat=x.interval)
    # no x.trans
    p <- ggOptScaleAxis(p, axis="x", scale=x.scale, breaks=x.breaks, verbose=verbose)
  }
  
  p <- ggOptCoordCartesian(p, df, x.id, y.id, x.lim.cart=x.lim.cart, y.lim.cart=y.lim.cart, 
                           coord.flip=coord.flip, verbose=verbose)
  
  p <- ggOptPalette(p, scale.to="fill", palette=palette)
  
  p <- ggOptLegend(p, legend.title=legend.title, legend.col=legend.col, legend.row=legend.row)
  
  p <- ggLabTitle(p, x.id, y.id, title=title, x.lab=x.lab, y.lab=y.lab)
  if (no.panel.border)
    p <- ggThemeAxis(p, title.size=title.size)
  else 
    p <- ggThemePanelBorder(p, title.size=title.size)
  
  p <- ggThemeOthers(p, x.text.angle=x.text.angle, legend.position=legend.position, 
                     legend.direction=legend.direction)
  
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
#' box.plot <- ggBoxWhiskersPlot(df, x.id="test", y.id="performance")
#' 
#' @rdname ggPlot
ggBoxWhiskersPlot <- function(df, x.id, y.id, fill.id=NULL, 
                              outlier.colour=alpha("black", 0.3), dodge.width=0.8,
                              x.facet.id=NULL, y.facet.id=NULL, coord.flip=FALSE,
                              y.trans="identity", auto.scale.y=FALSE, 
                              x.lim.cart=NULL, y.lim.cart=NULL, palette=NULL, 
                              legend.title=NULL, legend.col=1, legend.row=0, 
                              title="Box Whiskers Plot", title.size = 10, 
                              x.lab="x.id", y.lab="y.id", x.text.angle=0, 
                              no.panel.border=FALSE, verbose=TRUE) {
  p <- ggInit(df=df, x.id=x.id, y.id=y.id, fill.id=fill.id)
  if (! is.null(fill.id)) 
    p <- p + geom_boxplot(outlier.colour=outlier.colour, position=position_dodge(width=dodge.width))
  else 
    p <- p + geom_boxplot(outlier.colour=outlier.colour)
  
  p <- p + scale_shape(solid = FALSE) #+ geom_jitter(alpha = 0.5) 
  
  col.names <- colnames(df)
  p <- ggOptFacetGrid(p, col.names, x.facet.id=x.facet.id, y.facet.id=y.facet.id)
  
  if (auto.scale.y) {
    y.max <- max(df[,y.id])
    p <- ggOptScaleAxis(p, axis="y", scale="continuous", trans=y.trans, auto.scale.max=y.max)
  } else {
    p <- ggOptScaleAxis(p, axis="y", scale="continuous", trans=y.trans)
  }
  
  p <- ggOptCoordCartesian(p, df, x.id, y.id, x.lim.cart=x.lim.cart, y.lim.cart=y.lim.cart, 
                           coord.flip=coord.flip, verbose=verbose)
  
  p <- ggOptPalette(p, scale.to="fill", palette=palette)
  
  p <- ggOptLegend(p, legend.title=legend.title, legend.col=legend.col, legend.row=legend.row)
  
  p <- ggLabTitle(p, x.id, y.id, title=title, x.lab=x.lab, y.lab=y.lab)
  if (no.panel.border)
    p <- ggThemeAxis(p, title.size=title.size)
  else 
    p <- ggThemePanelBorder(p, title.size=title.size)
  
  p <- ggThemeOthers(p, x.text.angle=x.text.angle)
  
  return(p)
}


#' @details 
#' \code{gtScatterPlot} uses one-line function to plot many types of scatter chart.
#' 
#' @param point.size The size of points from \code{\link{point.size}}. Default to 3.
#' @param colour.id,shape.id,link.id The column name in \code{df} to 
#' define how the data points are coloured, shaped, or linked according their values.
#' @param ellipsed.id The column name in \code{df} to define 
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
#' pdf.gtable(g.table, fig.path="clusters-scatter-plot.pdf") 
#' 
#' @rdname ggPlot
gtScatterPlot <- function(df, x.id, y.id, colour.id=NULL, shape.id=NULL, 
                          shapes=NULL, point.size=3, x.facet.id=NULL, y.facet.id=NULL, 
                          link.id=NULL, ellipsed.id=NULL, text.id=NULL, 
                          text.data = NULL, text.size = 3, text.hjust=-0.1, 
                          text.vjust = -0.2, text.alpha = 0.5, coord.flip=FALSE,
                          xintercept=NULL, yintercept=NULL, line.type=2,
                          x.lim.cart=NULL, y.lim.cart=NULL, palette=NULL,
                          legend.title=NULL, legend.col=1, legend.row=0,  
                          title="Scatter Plot", title.size = 10, x.lab="x.id", y.lab="y.id", 
                          no.panel.border=FALSE, verbose=TRUE) {
  p <- ggInit(df=df, x.id=x.id, y.id=y.id, colour.id=colour.id)
  col.names <- colnames(df)
  
  if (! is.null(shape.id) && is.null(shapes)) {
    # The shape palette can deal with a maximum of 6 discrete values
    n_shape <- length(unique(df[,shape.id]))
    shapes <- seq(1, (1 + n_shape-1))
  }
  p <- ggOptPointAndShape(p, col.names, shape.id=shape.id, shapes=shapes, point.size=point.size)
  
  p <- ggOptFacetGrid(p, col.names, x.facet.id=x.facet.id, y.facet.id=y.facet.id)
  
  p <- ggOptEllipse(p, col.names, ellipsed.id=ellipsed.id)
  
  if (! is.null(link.id)) {
    suppressMessages(require(data.table))
    # Convex hull http://stackoverflow.com/questions/16428962/convex-hull-ggplot-using-data-tables-in-r
    df.dt <- data.table(df, key = link.id)
    chull.txt <- paste0('df.dt[, .SD[chull(x.id, y.id)], by = link.id ]')
    cat("chull.cmd : ", chull.txt, "\n")
    # hulls <- df.dt[, .SD[chull(MDS1, MDS2)], by = link.id]
    chull.cmd <- parse(text = chull.txt) 
    hulls <- eval(chull.cmd)
    
    p <- p + geom_polygon(data = hulls, aes_string(mapping=link.id), fill = NA, alpha = 0.5)
  }
  
  p <- ggOptText(p, col.names, text.id=text.id, text.data=text.data, colour.id=colour.id, 
                 text.size=text.size, text.hjust=text.hjust, text.vjust=text.vjust, 
                 text.alpha=text.alpha)
  
  if (! is.null(xintercept))
    p <- p + geom_vline(xintercept=xintercept,linetype=line.type)
  if (! is.null(yintercept))
    p <- p + geom_hline(yintercept=yintercept,linetype=line.type) 
  
  p <- ggOptCoordCartesian(p, df, x.id, y.id, x.lim.cart=x.lim.cart, y.lim.cart=y.lim.cart, 
                           coord.flip=coord.flip, verbose=verbose)
  
  p <- ggOptPalette(p, palette=palette)
  
  p <- ggOptLegend(p, legend.title=legend.title, legend.col=legend.col, legend.row=legend.row)
  
  p <- ggLabTitle(p, x.id, y.id, title=title, x.lab=x.lab, y.lab=y.lab)
  if (no.panel.border)
    p <- ggThemeAxis(p, title.size=title.size)
  else 
    p <- ggThemePanelBorder(p, title.size=title.size)
  
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
gtLine <- function(df, x.id, y.id, group.id=NULL, colour.id=NULL,  
                   line.size=0.5, line.type = 2, line.alpha=0.75, 
                   shape.id=NULL, shapes=NULL, point.size=3, point.data=NULL,
                   x.facet.id=NULL, y.facet.id=NULL, coord.flip=FALSE,
                   x.trans="identity", auto.scale.x=FALSE, y.trans="identity", auto.scale.y=FALSE,
                   text.id=NULL, text.data = NULL, text.size = 3, 
                   text.hjust=-0.1, text.vjust = -0.2, text.alpha = 0.5, 
                   x.lim.cart=NULL, y.lim.cart=NULL, palette=NULL, 
                   legend.title=NULL, legend.col=1, legend.row=0, x.text.angle=0, 
                   title="", title.size = 10, x.lab="x.id", y.lab="y.id", 
                   no.panel.border=FALSE, verbose=TRUE) {
  p <- ggInit(df=df, x.id=x.id, y.id=y.id, group.id=group.id, colour.id=colour.id)
  col.names <- colnames(df)
  
  p <- p + geom_line(size=line.size, linetype=line.type, alpha=line.alpha) 
  
  if (! is.null(shape.id) && is.null(shapes)) {
    # The shape palette can deal with a maximum of 6 discrete values
    n_shape <- length(unique(df[,shape.id]))
    shapes <- seq(1, (1 + n_shape-1))
  }
  p <- ggOptPointAndShape(p, col.names, shape.id=shape.id, data=point.data, 
                          shapes=shapes, point.size=point.size)
  
  if (auto.scale.x) {
    x.max <- max(df[,x.id])
    p <- ggOptScaleAxis(p, axis="x", scale="continuous", trans=x.trans, auto.scale.max=x.max)
  } else {
    p <- ggOptScaleAxis(p, axis="x", scale="continuous", trans=x.trans)
  }
  if (auto.scale.y) {
    y.max <- max(df[,y.id])
    p <- ggOptScaleAxis(p, axis="y", scale="continuous", trans=y.trans, auto.scale.max=y.max)
  } else {
    p <- ggOptScaleAxis(p, axis="y", scale="continuous", trans=y.trans)
  }
  
  p <- ggOptText(p, col.names, text.id=text.id, text.data=text.data, colour.id=colour.id, 
                 text.size=text.size, text.hjust=text.hjust, text.vjust=text.vjust, 
                 text.alpha=text.alpha)
  
  p <- ggOptCoordCartesian(p, df, x.id, y.id, x.lim.cart=x.lim.cart, y.lim.cart=y.lim.cart, 
                           coord.flip=coord.flip, verbose=verbose)
  
  p <- ggOptPalette(p, palette=palette)
  
  p <- ggOptLegend(p, legend.title=legend.title, legend.col=legend.col, legend.row=legend.row)
  
  p <- ggLabTitle(p, x.id, y.id, title=title, x.lab=x.lab, y.lab=y.lab)
  if (no.panel.border)
    p <- ggThemeAxis(p, title.size=title.size)
  else 
    p <- ggThemePanelBorder(p, title.size=title.size)
  
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  
  return(gt)
}

