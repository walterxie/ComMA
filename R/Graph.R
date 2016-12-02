# Graph
# Author: Walter Xie
# Accessed on 12 Apr 2016

#' @name ggPlot
#' @title One-line function to make \pkg{ggplot} 
#'
#' @description 
#' Simplify \pkg{ggplot} graphic coding into one-line functions that provides frequently 
#' used graphies in a convenient way. They are very basic charts here, and can be used or 
#' extended into various common charts, 
#' such as percentage bar chart \code{\link{ggPercentageBarChart}}, 
#' nonmetric multidimensional scaling plot \code{\link{gtNMDSPlot}}, 
#' PCA plot \code{\link{gtPCAPlot}}, 
#' rarefaction curves \code{\link{gtRarefactionCurve}}.
#' 
#' And also some useful charts that we derived from our publications, 
#' such as group abundance bar chart \code{\link{ggGroupAbundanceBar}},
#' and Y across X bar chart \code{\link{ggGroupAbundanceBar}}. 
#' Both refer to \url{http://dx.doi.org/10.1186/s13742-015-0086-1}.
#'
#' @return 
#' If the function returns a \code{\link{ggplot}} object, then its name starts with "gg". 
#' It needs to use \code{\link{pdf.ggplot}} to create pdf. 
#' It also keeps the expandability using '+'.  
#' 
#' If the function returns a \code{\link{gtable}} object, then its name starts with "gt".
#' This kind of functions use \code{\link{unclip.ggplot}} to turns off clipping for a 
#' \code{\link{ggplot}} object, but returns a \code{\link{gtable}} object.
#' It needs to use \code{\link{pdf.gtplot}} to create pdf. And also \code{\link{plot.gtable}} 
#' simplifies the code to plot gtable object in console.
#' @note 
#' All basic charts are designed to return a \code{\link{ggplot}} object for easy 
#' extension, you may need to turn off clipping after call, such as 
#' \code{ggScatterPlot} and \code{ggLineWithPoints}.
#' 
#' @param df A data frame used for plot. 
#' @param df.to.melt A data frame required to \code{\link{melt}} 
#' before making a \pkg{ggplot} object, such as input of \code{ggHeatmap}. 
#' At least one column should be \code{melt.id}, and the other columns should be values, 
#' unless they are fill.id, group.id, colour.id, etc. 
#' \strong{Note:} any extra column will disturb the result of \code{\link{melt}} function.
#'  
#' If row.names is going to be \code{melt.id}, 
#' then it can be inserted into the data frame before the plot. 
#' 
#' The data example 
#' \url{https://github.com/walterxie/ComMA/blob/master/data-raw/reads.phyla.txt}:
#' \tabular{rrrrr}{
#'   Phyla \tab 16s \tab 18s \tab ITS \tab TaxaGroup\cr
#'   Actinobacteria \tab 958 \tab 1 \tab 3 \tab Bacteria\cr
#'   Crenarchaeota \tab 1 \tab 0 \tab 0 \tab Archaea\cr
#'   Ascomycota \tab 2 \tab 765 \tab 971 \tab Fungi 
#' } 
#' @param melt.id A column name to \code{\link{melt}} 
#' and used as a \code{\link{factor}}, such as "plot" column.
#' @param x.id,y.id,fill.id,group.id The string of column names in \code{df} or 
#' \code{df.to.melt}, which use for \code{x, y, fill, group} in \code{\link{aes}}.
#' @param x.facet.id, y.facet.id The string of column names in \code{df},
#' which creates facets (a formula) in \code{\link{facet_grid}}.
#' @param x.trans,y.trans The string defines the data scale used in either x-axis or y-axis, 
#' which can be "identity" standing for normal, or "per" standing for percentage, 
#' moreover either the name of a transformation object for \code{\link{scale_x_continuous}}
#' or \code{\link{scale_y_continuous}} (e.g. \code{trans="log"}), or the object itself. 
#' Built-in transformations include "asn", "atanh", "boxcox", "exp", "identity", 
#' "log", "log10", "log1p", "log2", "logit", "probability", "probit", "reciprocal", 
#' "reverse" and "sqrt". Default to "identity". 
#' @param x.lim.cart,y.lim.cart Setting limits on the coordinate system will zoom the plot, 
#' and will not change the underlying data like setting limits on a scale will. 
#' Refer to \code{\link{coord_cartesian}}. 
#' Set lower bound only to y-axis using y.lim.cart=c(1000,NA). Default NULL. 
#' @param title Graph title, set title="" to remove it from the plot.
#' @param x.lab,y.lab The label of x-axis or y-axis, such as plot names. 
#' Set x.lab="" to remove x-axis label from the plot. Default to NULL to do nothing.
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
#' Set palette=c("red","darkgreen","darkgrey",""), if \code{\link{scale_fill_manual}} 
#' has to be used when the actual length(palette) <=3. 
#' Default NULL to use \code{\link{ggplot}} default colours.  
#' @param text.id Label the data points according \code{text.id} column, 
#' such as "Row.names" column after \code{\link{merge}}.
#' @param text.size,text.hjust,text.vjust,text.alpha 
#' The arguments to adjust text in \code{\link{geom_text}} in the line or scatter plot.
#' @param legend.title.fill,legend.title.colour,legend.title.shape,legend.title.group,
#' legend.title.size The title of legend created by fill, colour, shape, group, or size. 
#' Set legend.title.*="" to remove legend.
#' @param no.legend Turning off some legends, such as, fill, shape, colour. 
#' Default to NULL, which is to do nothing.
#' @param legend.col,legend.row Customize the number of columns or rows for legend in bar chart. 
#' They cannot be used at the same time. Default not to use them, legend.col=1, legend.row=0. 
#' @param x.text,y.text If FALSE, then hide x or y axis labels in plot. 
#' Default to TRUE.
#' @param no.panel.border Add panel border or not. Default to FALSE.


#' @details 
#' \code{ggBarChart} is an one-line function to plot many types of bar chart, 
#' such as normal bars, log-scaled bars, percentage bars, and also grouping.
#' Refer to \code{\link{geom_bar}}.
#' 
#' @param bar.pos Position adjustment for bars, either as a string, 
#' or the result of a call to a position adjustment function. 
#' Default to "dodge" in \code{ggBarChart}, "stack" in \code{ggHistogram}. 
#' Use \code{fill.id} to generate group percentage bars.
#' @param bar.stat Determine what is mapped to bar height in \code{ggBarChart}. 
#' Default to "identity", 
#' which defines the heights of the bars to represent values in the data.
#' Refer to \code{\link{geom_bar}}.
#' @keywords graph
#' @export
#' @examples
#' data(model.test)
#' model.test.mac <- model.test[model.test$OS=="Mac",]
#' ggBarChart(model.test.mac, x.id="test", y.id="performance", fill.id="model", x.text.angle=90)
#' ggBarChart(model.test.mac, x.id="test", y.id="performance", fill.id="model", x.text.angle=90, bar.pos="stack")
#' ggBarChart(model.test.mac, x.id="test", y.id="performance", fill.id="model", x.text.angle=90, bar.pos="fill", y.trans="per")
#' 
#' @rdname ggPlot
ggBarChart <- function(df, x.id, y.id, fill.id=NULL, 
                       bar.pos="dodge", bar.stat="identity", x.interval=0, 
                       x.trans="identity", x.scale="discrete", auto.scale.x=FALSE, 
                       y.trans="identity", y.scale="continuous", auto.scale.y=FALSE,
                       x.facet.id=NULL, y.facet.id=NULL, facet.scales="fixed", 
                       facet.space="fixed", facet.shrink=TRUE, facet.drop = TRUE,
                       xintercept=NULL, yintercept=NULL, line.type=2,
                       x.lim.cart=NULL, y.lim.cart=NULL, palette=NULL,
                       legend.title=NULL, legend.col=1, legend.row=0, 
                       title="Bar Chart", title.size=10, x.lab=NULL, y.lab=NULL, 
                       no.legend=NULL, legend.position="right", legend.direction="vertical",
                       coord.flip=FALSE, x.text.angle=0, x.text=TRUE, y.text=TRUE, 
                       plot.margin.cm=NULL, no.panel.border=FALSE, verbose=TRUE) {
  p <- ggInit(df=df, x.id=x.id, y.id=y.id, fill.id=fill.id, verbose=verbose)
  p <- p + geom_bar(position=bar.pos, stat=bar.stat) 
  
  col.names <- colnames(df)
  p <- ggOptFacetGrid(p, col.names, x.facet.id=x.facet.id, y.facet.id=y.facet.id, 
                      scales=facet.scales, space=facet.space, shrink=facet.shrink, 
                      drop=facet.drop, verbose=verbose)
  
  if (auto.scale.x) {
    x.max <- max(df[,x.id])
    p <- ggOptScaleAxis(p, axis="x", scale=x.scale, trans=x.trans, 
                        auto.scale.max=x.max, verbose=verbose)
  } else if (x.interval > 0) {
    #x.breaks <- seq(min(df[,x.id]), max(df[,x.id]), x.interval)
    x.breaks <- window(unique(df[,x.id]), deltat=x.interval)
    # no x.trans
    p <- ggOptScaleAxis(p, axis="x", scale=x.scale, breaks=x.breaks, verbose=verbose)
  } else {
    p <- ggOptScaleAxis(p, axis="x", scale=x.scale, trans=x.trans, 
                        verbose=verbose)
  }
  if (auto.scale.y) {
    y.max <- max(df[,y.id])
    p <- ggOptScaleAxis(p, axis="y", scale=y.scale, trans=y.trans, 
                        auto.scale.max=y.max, verbose=verbose)
  } else {
    p <- ggOptScaleAxis(p, axis="y", scale=y.scale, trans=y.trans, 
                        verbose=verbose)
  }
  
  p <- ggOptCoordCartesian(p, df, x.id, y.id, x.lim.cart=x.lim.cart, y.lim.cart=y.lim.cart, 
                           coord.flip=coord.flip, verbose=verbose)
  
  p <- ggOptPalette(p, scale.to="fill", palette=palette, verbose=verbose)
  
  p <- ggOptLegend(p, legend.title.fill=legend.title, no.legend=no.legend, 
                   legend.col=legend.col, legend.row=legend.row)
  
  p <- ggLabTitle(p, x.id, y.id, title=title, x.lab=x.lab, y.lab=y.lab)
  if (no.panel.border)
    p <- ggThemeAxis(p, title.size=title.size)
  else 
    p <- ggThemePanelBorder(p, title.size=title.size)
  
  p <- ggThemeOthers(p, x.text.angle=x.text.angle, legend.position=legend.position, 
                     legend.direction=legend.direction, x.text=x.text, y.text=y.text, 
                     plot.margin.cm=plot.margin.cm, verbose=verbose)
  
  return(p)
}

#' @details 
#' \code{ggHistogram} is an one-line function to plot a 1d distribution by dividing 
#' into bins and counting the number of observations in each bin.
#' Refer to \code{\link{geom_histogram}}.
#' 
#' @param bins Number of bins for \code{ggHistogram}. Overridden by \code{binwidth}. 
#' Defaults to 30.
#' or the result of a call to a position adjustment function. Default to "stack". 
#' Use \code{fill.id} to generate group percentage bars.
#' @param binwidth The width of the bins for \code{ggHistogram}. 
#' The default is to use \code{bins} that cover the range of the data. 
#' You should always override this value, exploring multiple widths to find 
#' the best to illustrate the stories in your data. 
#' The bin width of a date variable is the number of days in each time; 
#' the bin width of a time variable is the number of seconds.
#' Refer to \code{\link{geom_histogram}}.
#' @keywords graph
#' @export
#' @examples
#' ggHistogram(model.test.mac, x.id="performance", fill.id="model", x.text.angle=90)
#'  
#' @rdname ggPlot
ggHistogram <- function(df, x.id, fill.id=NULL, 
                        bar.pos="stack", bins = NULL, binwidth = NULL, x.interval=0,
                        x.trans="identity", x.scale="continuous", auto.scale.x=FALSE, 
                        y.trans="identity",
                        x.facet.id=NULL, y.facet.id=NULL, coord.flip=FALSE,
                        xintercept=NULL, yintercept=NULL, line.type=2,
                        x.lim.cart=NULL, y.lim.cart=NULL, palette=NULL,
                        legend.title=NULL, legend.col=1, legend.row=0, 
                        title="Histogram", title.size=10, x.lab=NULL, y.lab=NULL, 
                        no.legend=NULL, legend.position="right", legend.direction="vertical",
                        x.text.angle=0, x.text=TRUE, y.text=TRUE, 
                        plot.margin.cm=NULL, no.panel.border=FALSE, verbose=TRUE) {
  p <- ggInit(df=df, x.id=x.id, fill.id=fill.id, verbose=verbose)
  p <- p + geom_histogram(stat="bin", position=bar.pos, bins=bins, binwidth=binwidth) 
  
  col.names <- colnames(df)
  p <- ggOptFacetGrid(p, col.names, x.facet.id=x.facet.id, 
                      y.facet.id=y.facet.id, verbose=verbose)
  
  if (auto.scale.x) {
    x.max <- max(df[,x.id])
    p <- ggOptScaleAxis(p, axis="x", scale=x.scale, trans=x.trans, 
                        auto.scale.max=x.max, verbose=verbose)
  } else if (x.interval > 0) {
    #x.breaks <- seq(min(df[,x.id]), max(df[,x.id]), x.interval)
    x.breaks <- window(unique(df[,x.id]), deltat=x.interval)
    # no x.trans
    p <- ggOptScaleAxis(p, axis="x", scale=x.scale, breaks=x.breaks, verbose=verbose)
  } else {
    p <- ggOptScaleAxis(p, axis="x", scale=x.scale, trans=x.trans, 
                        verbose=verbose)
  }
  
  p <- ggOptScaleAxis(p, axis="y", scale="continuous", trans=y.trans, verbose=verbose)
 
  
  p <- ggOptCoordCartesian(p, df, x.id, y.id, x.lim.cart=x.lim.cart, y.lim.cart=y.lim.cart, 
                           coord.flip=coord.flip, verbose=verbose)
  
  p <- ggOptPalette(p, scale.to="fill", palette=palette, verbose=verbose)
  
  p <- ggOptLegend(p, legend.title.fill=legend.title, no.legend=no.legend, 
                   legend.col=legend.col, legend.row=legend.row)
  
  p <- ggLabTitle(p, x.id, y.id, title=title, x.lab=x.lab, y.lab=y.lab)
  if (no.panel.border)
    p <- ggThemeAxis(p, title.size=title.size)
  else 
    p <- ggThemePanelBorder(p, title.size=title.size)
  
  p <- ggThemeOthers(p, x.text.angle=x.text.angle, legend.position=legend.position, 
                     legend.direction=legend.direction, x.text=x.text, y.text=y.text, 
                     plot.margin.cm=plot.margin.cm, verbose=verbose)
  
  return(p)
}

#' @details 
#' \code{ggScatterPlot} uses one-line function to plot many types of scatter chart.
#' Refer to \code{\link{geom_point}}.
#' 
#' @param point.size,point.alpha The feature of points for \code{\link{geom_point}}. 
#' Default size to 3, alpha to 1.
#' @param text.or.point If 1, then display the text only; if 2, the default, then only points; 
#' if 3, then both the text and points.
#' @param colour.id,shape.id,link.id The column name in \code{df} to 
#' define how the data points are coloured, shaped, or linked according their values.
#' @param ellipsed.id The column name in \code{df} to define 
#' how to draw ellipse over data points, which is normally same as 
#' \code{colour.id} to show clusters.
#' @param shapes Manually define the shapes of points. Refer to \code{\link{scale_shape_manual}}, 
#' and \url{http://www.cookbook-r.com/Graphs/Shapes_and_line_types/}.
#' @param scale.limits.min Manually set the minimum data range of the scale 
#' given colours to \code{c(scale.limits.min, max(df[,colour.id]))}. 
#' Refer to \code{limits} in \code{\link{discrete_scale}}.
#' @param xintercept,yintercept,linetype Add horizontal or vertical line. 
#' Refer to \code{\link{geom_hline}} or \code{\link{geom_vline}}.
#' @param text.repel Apply \code{\link{geom_text_repel}} in \pkg{ggrepel}.
#' It is better than \code{text.avoid.overlap} which removes overlapped text.
#' Default to FALSE.
#' @param text.avoid.overlap If TRUE, text that overlaps previous text 
#' in the same layer will not be plotted. Use it carefully. Default to FALSE.
#' @keywords graph
#' @export
#' @examples 
#' # plot 2 clusters
#' df.clusters <- random2Clusters()
#' df.clusters$labels <- rownames(df.clusters)
#' df.clusters
#' gg.plot <- ggScatterPlot(df.clusters, x.id="x", y.id="y", colour.id="group", shape.id="group",   
#'                          xintercept=0, yintercept=0, title="Clusters", palette="Set1")
#' # turns off clipping
#' g.table <- unclip.ggplot(gg.plot) 
#' plot(g.table)
#' 
#' # selective labeling for points x > 3 and y > 6
#' ggScatterPlot(df.clusters, x.id="x", y.id="y", colour.id="group", ellipsed.id="group",
#'               text.id="labels", text.data=subset(df.clusters, x > 3 & y > 6), 
#'               xintercept=0, yintercept=0, title="Clusters", palette="Set1")
#'  
#' @rdname ggPlot
ggScatterPlot <- function(df, x.id, y.id, colour.id=NULL, text.colour.id=NULL, 
                          shape.id=NULL, shapes=NULL, link.id=NULL, ellipsed.id=NULL, 
                          text.id=NULL, text.or.point=2, point.size=3, point.alpha=1,
                          text.data = NULL, text.size = 3, text.alpha = 0.5, 
                          text.hjust=-0.1, text.vjust = -0.2, text.avoid.overlap = FALSE,
                          text.repel = FALSE, box.padding = unit(0.25, "lines"), 
                          point.padding = unit(1e-06, "lines"), arrow = NULL, force = 1,
                          x.trans="identity", x.scale="continuous", auto.scale.x=FALSE, 
                          y.trans="identity", y.scale="continuous", auto.scale.y=FALSE,
                          x.facet.id=NULL, y.facet.id=NULL, facet.scales="fixed", 
                          facet.space="fixed", facet.shrink=TRUE, facet.drop = TRUE,
                          xintercept=NULL, yintercept=NULL, line.type=2,
                          x.lim.cart=NULL, y.lim.cart=NULL, palette=NULL, scale.limits.min=NULL, 
                          legend.title.colour=NULL, legend.title.shape=NULL,
                          legend.title.size=NULL, legend.col=1, legend.row=0,  
                          title="Scatter Plot", title.size = 10, x.lab=NULL, y.lab=NULL, 
                          no.legend=NULL, legend.position="right", legend.direction="vertical",
                          coord.flip=FALSE, x.text.angle=0, x.text=TRUE, y.text=TRUE, 
                          no.panel.border=FALSE, verbose=TRUE) {
  p <- ggInit(df=df, x.id=x.id, y.id=y.id, colour.id=colour.id, verbose=verbose)
  col.names <- colnames(df)
  
  # if text.id not NULL, make default text.or.point inactived
  if (!is.null(text.id) && text.or.point==2)
    text.or.point=3
  
  if (text.or.point > 1) {
    if (! is.null(shape.id) && is.null(shapes)) {
      # The shape palette can deal with a maximum of 6 discrete values
      n_shape <- length(unique(df[,shape.id]))
      shapes <- seq(1, (1 + n_shape-1))
    }
    p <- ggOptPointAndShape(p, col.names, shape.id=shape.id, shapes=shapes, 
                            point.size=point.size, point.alpha=point.alpha)
  }
  
  if (text.or.point != 2) {
    if (is.null(text.id)) 
      stop("Please specify 'text.id' to display text !")
    if (text.or.point == 1) {
      # prevent Error: Aesthetics must be either length 1 or the same as the data (5): size
      if (is.null(text.data))
        text.data=df
      if (text.hjust == -0.1)
        text.hjust=0 
      if (text.vjust == -0.2)
        text.vjust=0
      if (text.alpha == 0.5)
        text.alpha = 1
    }
    if (is.null(text.colour.id))
      text.colour.id <- colour.id
    p <- ggOptText(p, col.names, text.id=text.id, text.data=text.data, colour.id=text.colour.id, 
                   text.repel=text.repel, text.size=text.size, text.alpha=text.alpha, 
                   text.hjust=text.hjust, text.vjust=text.vjust, text.avoid.overlap=text.avoid.overlap, 
                   box.padding=box.padding, point.padding=point.padding, arrow=arrow, 
                   force=force, verbose=verbose)
  }
  
  p <- ggOptFacetGrid(p, col.names, x.facet.id=x.facet.id, y.facet.id=y.facet.id, 
                      scales=facet.scales, space=facet.space, shrink=facet.shrink, 
                      drop=facet.drop, verbose=verbose)
  
  p <- ggOptEllipse(p, col.names, ellipsed.id=ellipsed.id)
  
  if (! is.null(link.id)) {
    # TODO: not working properly
    # suppressMessages(require(data.table)) 
    # # Convex hull http://stackoverflow.com/questions/16428962/convex-hull-ggplot-using-data-tables-in-r
    # df.dt <- data.table(df, key = link.id)
    # chull.txt <- paste0('df.dt[, .SD[chull(x.id, y.id)], by = link.id ]')
    # cat("chull.cmd : ", chull.txt, "\n")
    # # hulls <- df.dt[, .SD[chull(MDS1, MDS2)], by = link.id]
    # chull.cmd <- parse(text = chull.txt) 
    # hulls <- eval(chull.cmd)
    # p <- p + geom_polygon(data = hulls, aes_string(mapping=link.id), fill = NA, alpha = 0.5)
    
    p <- p + geom_polygon(data = df, aes_string(mapping=link.id), fill = NA, alpha = 0.75)
  }
  
  if (auto.scale.x) {
    x.max <- max(df[,x.id])
    p <- ggOptScaleAxis(p, axis="x", scale=x.scale, trans=x.trans, 
                        auto.scale.max=x.max, verbose=verbose)
  } else {
    p <- ggOptScaleAxis(p, axis="x", scale=x.scale, trans=x.trans, 
                        verbose=verbose)
  }
  if (auto.scale.y) {
    y.max <- max(df[,y.id])
    p <- ggOptScaleAxis(p, axis="y", scale=y.scale, trans=y.trans, 
                        auto.scale.max=y.max, verbose=verbose)
  } else {
    p <- ggOptScaleAxis(p, axis="y", scale=y.scale, trans=y.trans, 
                        verbose=verbose)
  }
  
  if (! is.null(xintercept))
    p <- p + geom_vline(xintercept=xintercept,linetype=line.type)
  if (! is.null(yintercept))
    p <- p + geom_hline(yintercept=yintercept,linetype=line.type) 
  
  p <- ggOptCoordCartesian(p, df, x.id, y.id, x.lim.cart=x.lim.cart, y.lim.cart=y.lim.cart, 
                           coord.flip=coord.flip, verbose=verbose)
  
  colour.limits <- NULL
  if (!is.null(scale.limits.min))
    colour.limits <- c(scale.limits.min, max(df[,colour.id]))
  p <- ggOptPalette(p, palette=palette, limits=colour.limits, verbose=verbose)
  
  p <- ggOptLegend(p, legend.title.colour=legend.title.colour,  
                   legend.title.shape=legend.title.shape, legend.title.size=legend.title.size,
                   no.legend=no.legend, legend.col=legend.col, legend.row=legend.row)
  
  p <- ggLabTitle(p, x.id, y.id, title=title, x.lab=x.lab, y.lab=y.lab)
  if (no.panel.border)
    p <- ggThemeAxis(p, title.size=title.size)
  else 
    p <- ggThemePanelBorder(p, title.size=title.size)
  
  p <- ggThemeOthers(p, x.text.angle=x.text.angle, legend.position=legend.position, 
                     legend.direction=legend.direction, x.text=x.text, y.text=y.text, 
                     verbose=verbose)
  return(p)
}

#' @details 
#' \code{ggLineWithPoints} uses one-line function to plot a line or group of lines.
#' Refer to \code{\link{geom_line}}.
#' 
#' @param line.type,line.size,line.alpha The feature of lines for \code{\link{geom_line}}.
#' Default line.type to 1, size to 0.5, alpha to 1.
#' @param line.or.point If 1, then display the line only; if 2, then only points; 
#' if 3, the default, then both the line and points.
#' @param point.data An option to add extra points in the plot, such as \code{gtLine}. 
#' The data to be displayed as points in \code{\link{geom_point}}, which can be 
#' the subset of points to draw the line. 
#' If NULL, the default, the data is inherited from the plot data as specified in the call.
#' @keywords graph
#' @export
#' @examples 
#' mcmc.log <- readMCMCLog("data-raw/star.beast.log")
#' mcmc.log$state <- as.double(rownames(mcmc.log))
#' names(mcmc.log)
#' gg.plot <- ggLineWithPoints(mcmc.log[,c("TreeHeight.Species", "state")], x.id="state", y.id="TreeHeight.Species")
#' # turns off clipping
#' g.table <- unclip.ggplot(gg.plot) 
#' plot(g.table)
#' 
#' ggLineWithPoints(mcmc.log[,c("TreeHeight.Species", "state")], x.id="state", y.id="TreeHeight.Species", line.or.point=1)
#' ggLineWithPoints(mcmc.log[,c("TreeHeight.Species", "state")], x.id="state", y.id="TreeHeight.Species", line.or.point=2, point.size=1)
#' 
#' @rdname ggPlot
ggLineWithPoints <- function(df, x.id, y.id, group.id=NULL, colour.id=NULL, 
                             shape.id=NULL, shapes=NULL, text.id=NULL,  
                             line.or.point=3, line.size=0.5, line.type = 1, line.alpha=1, 
                             text.or.point=2, point.size=3, point.alpha=1, point.data=NULL, 
                             text.data = NULL, text.size = 3, text.alpha = 0.5, 
                             text.hjust=-0.1, text.vjust = -0.2, text.avoid.overlap = FALSE,
                             text.repel = FALSE, box.padding = unit(0.25, "lines"), 
                             point.padding = unit(1e-06, "lines"), arrow = NULL, force = 1,
                             x.facet.id=NULL, y.facet.id=NULL, facet.scales="fixed", 
                             facet.space="fixed", facet.shrink=TRUE, facet.drop = TRUE,
                             x.trans="identity", x.scale="continuous", auto.scale.x=FALSE, 
                             y.trans="identity", y.scale="continuous", auto.scale.y=FALSE,
                             x.lim.cart=NULL, y.lim.cart=NULL, palette=NULL, 
                             legend.title.group=NULL, legend.title.colour=NULL, 
                             legend.title.shape=NULL, legend.col=1, legend.row=0, 
                             title="", title.size = 10, x.lab=NULL, y.lab=NULL, 
                             no.legend=NULL, legend.position="right", legend.direction="vertical",
                             coord.flip=FALSE, x.text.angle=0, x.text=TRUE, y.text=TRUE, 
                             no.panel.border=FALSE, verbose=TRUE) {
  p <- ggInit(df=df, x.id=x.id, y.id=y.id, group.id=group.id, colour.id=colour.id, verbose=verbose)
  col.names <- colnames(df)
  
  if (line.or.point != 2)
    p <- p + geom_line(size=line.size, linetype=line.type, alpha=line.alpha) 
  
  if (line.or.point > 1) {
    if (! is.null(shape.id) && is.null(shapes)) {
      # The shape palette can deal with a maximum of 6 discrete values
      n_shape <- length(unique(df[,shape.id]))
      shapes <- seq(1, (1 + n_shape-1))
    }
    p <- ggOptPointAndShape(p, col.names, shape.id=shape.id, data=point.data, 
                            shapes=shapes, point.size=point.size, point.alpha=point.alpha)
  }
  
  if (auto.scale.x) {
    x.max <- max(df[,x.id])
    p <- ggOptScaleAxis(p, axis="x", scale=x.scale, trans=x.trans, 
                        auto.scale.max=x.max, verbose=verbose)
  } else {
    p <- ggOptScaleAxis(p, axis="x", scale=x.scale, trans=x.trans, 
                        verbose=verbose)
  }
  if (auto.scale.y) {
    y.max <- max(df[,y.id])
    p <- ggOptScaleAxis(p, axis="y", scale=y.scale, trans=y.trans, 
                        auto.scale.max=y.max, verbose=verbose)
  } else {
    p <- ggOptScaleAxis(p, axis="y", scale=y.scale, trans=y.trans, 
                        verbose=verbose)
  }
  
  p <- ggOptFacetGrid(p, col.names, x.facet.id=x.facet.id, y.facet.id=y.facet.id, 
                      scales=facet.scales, space=facet.space, shrink=facet.shrink, 
                      drop=facet.drop, verbose=verbose)
  
  p <- ggOptText(p, col.names, text.id=text.id, text.data=text.data, colour.id=text.colour.id, 
                 text.repel=text.repel, text.size=text.size, text.alpha=text.alpha, 
                 text.hjust=text.hjust, text.vjust=text.vjust, text.avoid.overlap=text.avoid.overlap, 
                 box.padding=box.padding, point.padding=point.padding, arrow=arrow, 
                 force=force, verbose=verbose)
  
  p <- ggOptCoordCartesian(p, df, x.id, y.id, x.lim.cart=x.lim.cart, y.lim.cart=y.lim.cart, 
                           coord.flip=coord.flip, verbose=verbose)
  
  p <- ggOptPalette(p, palette=palette, verbose=verbose)
  
  p <- ggOptLegend(p, legend.title.group=legend.title.group, 
                   legend.title.colour=legend.title.colour, legend.title.shape=legend.title.shape, 
                   no.legend=no.legend, legend.col=legend.col, legend.row=legend.row)
  
  p <- ggLabTitle(p, x.id, y.id, title=title, x.lab=x.lab, y.lab=y.lab)
  if (no.panel.border)
    p <- ggThemeAxis(p, title.size=title.size)
  else 
    p <- ggThemePanelBorder(p, title.size=title.size)
  
  p <- ggThemeOthers(p, x.text.angle=x.text.angle, legend.position=legend.position, 
                     legend.direction=legend.direction, x.text=x.text, y.text=y.text, 
                     verbose=verbose)
  return(p)
}

#' @details 
#' \code{ggHeatmap} creates a heat map using ggplot. It is devried and improved 
#' from the code  
#' \url{https://learnr.wordpress.com/2010/01/26/ggplot2-quick-heatmap-plotting/}. 
#' 
#' @param low,high,mid,na.value,midpoint,colours,values The colour range of heatmap.
#' Default to low="white", high="steelblue". 
#' Refer to \code{\link{scale_fill_gradient}} or \code{\link{scale_fill_gradient2}}
#' (if midpoint is given). 
#' @param log.scale.colour If TRUE, then use log scale to the colour of heat map.
#' Default to FALSE.
#' @param limits,breaks Refer to \code{\link{discrete_scale}}.
#' @param auto.breaks.length,breaks.digits Heatmap ledgend breaks automatically.
#' If \code{auto.breaks.length <= 1}, then use \code{breaks}, the default is 5.
#' If values are between 0 and 1, then make sure \code{breaks.digits} is set properly,  
#' such as 2 in this case. 
#' @param add.label,label.digits \code{add.label} is logical, default to TRUE, 
#' to add value into the tiles of heatmap. 
#' \code{label.digits} determines the digit of the label.
#' @param x.levels,y.levels The levels to order x or y.
#' @keywords graph
#' @export
#' @examples 
#' ranks.by.group <- data.frame(plot=c("Plot03","Plot02","Plot01"), `16s`=c(3,2,1), 
#'                              `18s`=c(1,2,3), ITS=c(2,1,3), check.names = F)
#' ranks.by.group
#' gg.plot <- ggHeatmap(ranks.by.group, melt.id="plot")
#' pdf.ggplot(gg.plot, fig.path="plot-prior-example-heatmap.pdf") 
#' 
#' ggHeatmap(corr.tri, melt.id="gene", title="", legend.title="Correlations", label.digits=2, breaks.digits=2)
#' 
#' @rdname ggPlot
ggHeatmap <- function(df.to.melt, melt.id, low="white", high="steelblue", mid = "white",
                      midpoint = NULL, colours=NULL, values=NULL, na.value="transparent",
                      limits = NULL, breaks = waiver(), add.label=T, label.digits=1, 
                      auto.breaks.length = 5, breaks.digits = 0,
                      title="Heatmap", title.size = 10, x.lab="", y.lab="",
                      log.scale.colour=FALSE, legend.title="Counts",
                      x.lim.cart=NULL, y.lim.cart=NULL, coord.flip=FALSE,
                      no.legend=NULL, legend.position="right", legend.direction="vertical",
                      x.text.angle=90, x.text=TRUE, y.text=TRUE,
                      x.levels=c(), y.levels=c(), 
                      no.panel.border=FALSE, verbose=TRUE) {
  if (!is.element(tolower(melt.id), tolower(colnames(df.to.melt))))
    stop("Data frame column names do NOT have \"", melt.id, "\" for melt function !")
  
  suppressMessages(require(reshape2))
  df.melt <- melt(df.to.melt, id=c(melt.id))
  if(anyNA(df.melt$value)) # rm all NA
    df.melt <- na.omit(df.melt)
  
  if (length(x.levels)>1) {
    if (length(x.levels) != length(unique(df.melt$variable)))
      warning("x.levels length != x unique values !")
    df.melt$variable <- factor(df.melt$variable, ordered = TRUE, levels = x.levels)
  }
  if (length(y.levels)>1) {
    if (length(y.levels) != length(unique(df.melt[,melt.id])))
      warning("y.levels length != y unique values !")
    df.melt[,melt.id] <- factor(df.melt[,melt.id], ordered = TRUE, levels = y.levels)
  } else {
    df.melt[,melt.id] <- factor(df.melt[,melt.id], levels=sort(unique(df.melt[,melt.id]), decreasing = T))
  }
  
  suppressMessages(require(ggplot2))
  # variable is all group names, such as "16S" or "FUNGI"
  # value is ranks for each group
  p <- ggplot(df.melt, aes_string(x="variable", y=melt.id)) + geom_tile(aes(fill=value)) 
  
  if (length(breaks) == 0 && auto.breaks.length > 1) {
    min.v <- min(df.melt$value)
    max.v <- max(df.melt$value)
    breaks <- round(seq(min.v, max.v, length.out = auto.breaks.length), digits = breaks.digits)
    
    if (log.scale.colour) {
      if (min.v == 0)
        min.log <- 0
      else 
        min.log <- log(min.v)
      if (max.v == 0)
        max.log <- 0
      else 
        max.log <- log(max.v)
      
      breaks <- round(exp(seq(min.log, max.log, length.out = auto.breaks.length)), 
                      digits = breaks.digits)
    }
  }
  if (log.scale.colour) {
    p <- p + scale_fill_gradient(trans='log', na.value=na.value, low=low, high=high, 
                                 name=legend.title, breaks=breaks, limits=limits) 
  } else {
    if (!is.null(midpoint)) {
      p <- p + scale_fill_gradient2(midpoint=midpoint, na.value=na.value, low=low, high=high, 
                                    mid =mid, name=legend.title, breaks=breaks, limits=limits)
    } else if (!is.null(colours)) {
      if (auto.breaks.length > 1 && !is.null(values) && length(values) < 1) {
        require(scales)
        if (min.v < 0 && max.v > 0) {
          values=rescale(c(min.v,0-.Machine$double.eps,0,0+.Machine$double.eps,max.v))
          if (! 0 %in% breaks) {
            breaks <- sort(c(0, breaks), decreasing = T)
          }
        } else {
          values=rescale(c(min.v,max.v))
        }
      }
      p <- p + scale_fill_gradientn(colours=colours, na.value=na.value, values=values, 
                                    name=legend.title, breaks=breaks, limits=limits)
    } else {
      p <- p + scale_fill_gradient(na.value=na.value, low=low, high=high, 
                                   name=legend.title, breaks=breaks, limits=limits) 
    }
  }
  cat("The final breaks are : ", paste(breaks, collapse = ","), ".\n")
  if (length(breaks) > 2 && length(unique(breaks)) <= 2)
    warning("breaks are not set properly ! Try breaks.digits = 2")
  
  if (add.label) {
    p <- p + geom_text(aes(label=paste(round(value, digits = label.digits))))
  }
  
  p <- ggOptCoordCartesian(p, df.melt, x.id="variable", y.id=melt.id, 
                           x.lim.cart=x.lim.cart, y.lim.cart=y.lim.cart, 
                           coord.flip=coord.flip, verbose=verbose)
  
  p <- ggLabTitle(p, "", "", title=title, x.lab=x.lab, y.lab=y.lab)
  if (no.panel.border)
    p <- ggThemeAxis(p, title.size=title.size)
  else 
    p <- ggThemePanelBorder(p, title.size=title.size)
  
  p <- ggThemeOthers(p, x.text.angle=x.text.angle, legend.position=legend.position, 
                     legend.direction=legend.direction, x.text=x.text, y.text=y.text, 
                     verbose=verbose)
  
  return(p) 
}

#' @details 
#' \code{ggBoxWhiskersPlot} creates box Whiskers plot. 
#' Refer to \code{\link{geom_boxplot}}.
#' 
#' @param outlier.colour The colour of outliers in box whiskers plot 
#' used for \code{outlier.colour} in \code{\link{geom_boxplot}}. 
#' Default to alpha("black", 0.3).
#' @param dodge.width Dodging width, when different to the width of 
#' the individual elements in box Whiskers plot. 
#' This is useful when you want to align narrow geoms with wider geoms. 
#' Refer to \code{\link{position_dodge}}.
#' @param lwd change the line thinkness. 
#' @param fatten make the median line thinner relative to the other lines.
#' @keywords graph
#' @export
#' @examples
#' ggBoxWhiskersPlot(df, x.id="test", y.id="performance", colour.id="OS", x.text.angle=90)
#' 
#' @rdname ggPlot
ggBoxWhiskersPlot <- function(df, x.id, y.id, fill.id=NULL, colour.id=NULL,
                              outlier.colour=alpha("black", 0.3), dodge.width=0.8,
                              x.facet.id=NULL, y.facet.id=NULL, facet.scales="fixed", 
                              facet.space="fixed", facet.shrink=TRUE, facet.drop = TRUE,
                              y.trans="identity", auto.scale.y=FALSE, 
                              x.lim.cart=NULL, y.lim.cart=NULL, lwd=1, fatten=1,
                              palette=NULL, scale.to="fill", scale.type=NULL,
                              legend.title.fill=NULL, legend.title.colour=NULL, 
                              no.legend=NULL, legend.col=1, legend.row=0, 
                              legend.position="right", legend.direction="vertical",
                              title="Box Whiskers Plot", title.size = 10, 
                              x.lab=NULL, y.lab=NULL, 
                              coord.flip=FALSE, x.text.angle=0, x.text=TRUE, y.text=TRUE,
                              no.panel.border=FALSE, verbose=TRUE) {
  p <- ggInit(df=df, x.id=x.id, y.id=y.id, fill.id=fill.id, colour.id=colour.id, verbose=verbose)
  if (! is.null(fill.id)) 
    p <- p + geom_boxplot(outlier.colour=outlier.colour, position=position_dodge(width=dodge.width),
                          lwd=lwd, fatten=fatten)
  else 
    p <- p + geom_boxplot(outlier.colour=outlier.colour, lwd=lwd, fatten=fatten)
  
  p <- p + scale_shape(solid = FALSE) #+ geom_jitter(alpha = 0.5) 
  
  col.names <- colnames(df)
  p <- ggOptFacetGrid(p, col.names, x.facet.id=x.facet.id, y.facet.id=y.facet.id, 
                      scales=facet.scales, space=facet.space, shrink=facet.shrink, 
                      drop=facet.drop, verbose=verbose)
  
  if (auto.scale.y) {
    y.max <- max(df[,y.id])
    p <- ggOptScaleAxis(p, axis="y", scale="continuous", trans=y.trans, 
                        auto.scale.max=y.max, verbose=verbose)
  } else {
    p <- ggOptScaleAxis(p, axis="y", scale="continuous", trans=y.trans, 
                        verbose=verbose)
  }
  
  p <- ggOptCoordCartesian(p, df, x.id, y.id, x.lim.cart=x.lim.cart, y.lim.cart=y.lim.cart, 
                           coord.flip=coord.flip, verbose=verbose)
  
  p <- ggOptPalette(p, scale.to=scale.to, scale.type=scale.type, palette=palette, verbose=verbose)
  
  p <- ggOptLegend(p, legend.title.colour=legend.title.colour, 
                   legend.title.fill=legend.title.fill, no.legend=no.legend, 
                   legend.col=legend.col, legend.row=legend.row)
  
  p <- ggLabTitle(p, x.id, y.id, title=title, x.lab=x.lab, y.lab=y.lab)
  if (no.panel.border)
    p <- ggThemeAxis(p, title.size=title.size)
  else 
    p <- ggThemePanelBorder(p, title.size=title.size)
  
  p <- ggThemeOthers(p, x.text.angle=x.text.angle, legend.position=legend.position, 
                     legend.direction=legend.direction, x.text=x.text, y.text=y.text, 
                     verbose=verbose)
  
  return(p)
}

#' @details 
#' \code{ggDensityEstimate} is an one-line function to plot kernel density estimate (KDE), 
#' useful for display the distribution of variables with underlying smoothness. 
#' Give a column name \code{fill.id} to colour the filled area, and \code{colour.id} 
#' to colour the curve. Use \code{density.pos} to produce different types of density. 
#' Use \code{y.id="..count.."} to produce a conditional density estimate, 
#' which counts (density * n) variable instead of the default density.
#' Refer to \code{\link{geom_density}}.
#' 
#' @param density.pos Position adjustment for \code{ggDensityEstimate}, either as a string, 
#' or the result of a call to a position adjustment function. Default to "identity". 
#' Use "stack" to produce stacked density plots, 
#' and "fill" to plot density estimate in percentage scale.
#' @param density.alpha Modify colour transparency when \code{fill.id} is assigned to 
#' \code{ggDensityEstimate}. Default to 0.1. Refer to \code{\link{alpha}}. 
#' @keywords graph
#' @export
#' @examples
#' # prepare log
#' b.log <- ComMA::readFile("data-raw/star.beast.log", row.names=NULL)
#' df.melt <- melt(b.log, id="Sample")
#' df.TreeHeight <- df.melt[grep("TreeHeight", df.melt[,"variable"]),]
#' 
#' ggDensityEstimate(df.TreeHeight, x.id="value", colour.id="variable")
#' ggDensityEstimate(df.TreeHeight, x.id="value", fill.id="variable", colour.id="variable")
#' # stacked density plot to lose marginal densities
#' ggDensityEstimate(df.TreeHeight, x.id="value", fill.id="variable", density.pos="stack", density.alpha=1)
#' # conditional density plot to preserve marginal densities
#' ggDensityEstimate(df.TreeHeight, x.id="value", y.id="..count..", fill.id="variable", density.pos="stack", density.alpha=1)
#' # percentage scale
#' ggDensityEstimate(df.TreeHeight, x.id="value", fill.id="variable", density.pos="fill", density.alpha=1)
#' 
#' @rdname ggPlot
ggDensityEstimate <- function(df, x.id, y.id=NULL, fill.id=NULL, colour.id=NULL,
                              density.pos="identity", density.alpha=0.1,
                              x.facet.id=NULL, y.facet.id=NULL, 
                              x.lim.cart=NULL, y.lim.cart=NULL,   
                              x.trans="identity", x.scale="continuous", auto.scale.x=FALSE, 
                              y.trans="identity", y.scale="continuous", auto.scale.y=FALSE,
                              fill.palette=NULL, colour.palette=NULL, coord.flip=FALSE,
                              legend.title.colour=NULL, legend.title.fill=NULL,
                              legend.col=1, legend.row=0, 
                              title="Kernel Density Estimate", title.size=10, 
                              x.lab=NULL, y.lab=NULL, no.legend=NULL, 
                              legend.position="right", legend.direction="vertical",
                              x.text.angle=0, x.text=TRUE, y.text=TRUE,
                              no.panel.border=FALSE, verbose=TRUE) {
  p <- ggInit(df=df, x.id=x.id, y.id=y.id, fill.id=fill.id, colour.id=colour.id, verbose=verbose)
  p <- p + geom_density(position=density.pos, alpha=density.alpha) 
  
  col.names <- colnames(df)
  p <- ggOptFacetGrid(p, col.names, x.facet.id=x.facet.id, 
                      y.facet.id=y.facet.id, verbose=verbose)
  
  if (auto.scale.x) {
    x.max <- max(df[,x.id])
    p <- ggOptScaleAxis(p, axis="x", scale=x.scale, trans=x.trans, 
                        auto.scale.max=x.max, verbose=verbose)
  } else {
    p <- ggOptScaleAxis(p, axis="x", scale=x.scale, trans=x.trans, 
                        verbose=verbose)
  }
  if (auto.scale.y) {
    y.max <- max(df[,y.id])
    p <- ggOptScaleAxis(p, axis="y", scale=y.scale, trans=y.trans, 
                        auto.scale.max=y.max, verbose=verbose)
  } else {
    p <- ggOptScaleAxis(p, axis="y", scale=y.scale, trans=y.trans, 
                        verbose=verbose)
  }
  
  p <- ggOptCoordCartesian(p, df, x.id, y.id, x.lim.cart=x.lim.cart, y.lim.cart=y.lim.cart, 
                           coord.flip=coord.flip, verbose=verbose)
  
  p <- ggOptPalette(p, scale.to="fill", palette=fill.palette, verbose=verbose)
  p <- ggOptPalette(p, scale.to="colour", palette=colour.palette, verbose=verbose)
  
  p <- ggOptLegend(p, legend.title.colour=legend.title.colour, 
                   legend.title.fill=legend.title.fill, no.legend=no.legend, 
                   legend.col=legend.col, legend.row=legend.row)
  
  if (density.pos=="stack" && title=="Kernel Density Estimate")
    title <- paste("Conditional", title)
  p <- ggLabTitle(p, x.id, y.id, title=title, x.lab=x.lab, y.lab=y.lab)
  if (no.panel.border)
    p <- ggThemeAxis(p, title.size=title.size)
  else 
    p <- ggThemePanelBorder(p, title.size=title.size)
  
  p <- ggThemeOthers(p, x.text.angle=x.text.angle, legend.position=legend.position, 
                     legend.direction=legend.direction, x.text=x.text, y.text=y.text, 
                     verbose=verbose)
  
  return(p)
}


