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
    stop(paste0("Data frame column names do NOT have \"", melt.id, "\" for melt function !"))
  
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
#' @param x.scale,y.scale The string defines the data scale used in either x-axis or y-axis, 
#' which can be "identity" standing for normal, or "per" standing for percentage, 
#' moreover either the name of a transformation object for \code{\link{scale_x_continuous}}
#' or \code{\link{scale_y_continuous}} (e.g. \code{trans="log"}), or the object itself. 
#' Built-in transformations include "asn", "atanh", "boxcox", "exp", "identity", 
#' "log", "log10", "log1p", "log2", "logit", "probability", "probit", "reciprocal", 
#' "reverse" and "sqrt". Default to "identity". 
#' @param x.lim.cart,y.lim.cart Setting limits on the coordinate system will zoom the plot, 
#' and will not change the underlying data like setting limits on a scale will. 
#' Refer to \code{\link{coord_cartesian}}. Set lower bound only to y-axis using y.lim.cart=c(1000,NA). 
#' Default missing. 
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
#' Default missing to use \code{\link{ggplot}} default colours.  
#' @param legend.col,legend.nrow Customize the number of columns or rows for legend in bar chart. 
#' They cannot be used at the same time. Default not to use them, legend.col=1, legend.nrow=0. 
#' @keywords graph
#' @export
#' @examples
#' # log-scale y
#' bar.chart <- ggBarChart(df.melt, x.id="test", y.id="seconds", fill.id="version", y.scale="log")
#' # percentage bars without grouping in one bar each
#' bar.chart <- ggBarChart(df.melt, x.id="test", y.id="percentage", fill.id="model", y.scale="per")
#' # percentage bars one group in one bar
#' bar.chart <- ggBarChart(df.melt, x.id="test", y.id="percentage", fill.id="model", bar.pos="fill", y.scale="per")
#' 
#' # the number of OTUs (y-axis) across the number of samples (x-axis)
#' communityMatrix <- getCommunityMatrix("16S", isPlot=TRUE, minAbund=1)
#' cm.aggre <- cmYAcrossX(communityMatrix)
#' df.melt <- melt(cm.aggre, id="samples")
#' bar.chart <- ggBarChart(df.melt, x.id="samples", y.id="value", fill.id="variable", y.scale="log", 
#'                         y.lab="", legend.title="", x.interval=1)
#' 
#' @rdname ggPlot
ggBarChart <- function(df.melt, x.id, y.id, fill.id, x.facet.id, y.facet.id, 
                       bar.pos="dodge", y.scale="identity", 
                       rotate.x.text=FALSE, x.interval=0, x.lim.cart, y.lim.cart,
                       palette, legend.col=1, legend.nrow=0, legend.title, 
                       title="Bar Chart", title.size = 10, x.lab="x.id", y.lab="y.id") {
  if (!is.element(tolower(x.id), tolower(colnames(df.melt))))
    stop(paste0("Data frame do NOT have column name \"", x.id, "\" !"))
  if (!is.element(tolower(y.id), tolower(colnames(df.melt))))
    stop(paste0("Data frame do NOT have column name \"", y.id, "\" !"))
  
  require(ggplot2)
  if (missing(fill.id)) {
    p <- ggplot(df.melt, aes_string(x = x.id, y = y.id)) 
  } else {
    if (!is.element(tolower(fill.id), tolower(colnames(df.melt))))
      stop(paste0("Data frame do NOT have column name \"", fill.id, "\" !"))
    
    p <- ggplot(df.melt, aes_string(x = x.id, y = y.id, fill = fill.id))
  }
  p <- p + geom_bar(position = bar.pos, stat = "identity") 
  
  if (! missing(x.facet.id) || ! missing(y.facet.id)) {
    fac1 <- "."
    fac2 <- "."
    if (! missing(x.facet.id))
      fac1 <- x.facet.id
    if (! missing(y.facet.id))
      fac2 <- y.facet.id
    
    p <- p + facet_grid(reformulate(fac2,fac1))
  }
  
  if (y.scale=="identity") {
    p <- p + scale_y_continuous(expand = c(0,0))
  } else if (y.scale=="per") {
    require(scales)
    p <- p + scale_y_continuous(expand = c(0,0), labels = percent_format()) 
  } else {
    y.breaks <- ComMA::get_breaks_positive_values(max(df.melt[,y.id], start=c(0)))
    p <- p + scale_y_continuous(trans=y.scale, expand = c(0,0), 
                                labels = ComMA::scientific_10, breaks = y.breaks) 
  }
  
  if (x.interval > 0) {
    x.breaks <- seq(min(df.melt[,x.id]), max(df.melt[,x.id]), x.interval)
    p <- p + scale_x_discrete(breaks=x.breaks) 
  }
  
  if (! missing(x.lim.cart)) {
    if (length(x.lim.cart) != 2)
      stop("Invalid format, use x.lim.cart = c(1000,NA) !")
    if (which(is.na(x.lim.cart))==1) {
      x.lim.cart[1] = min(df.melt[,x.id])
    } else if (which(is.na(x.lim.cart))==2) {
      x.lim.cart[2] = max(df.melt[,x.id])
    }
    p <- p + coord_cartesian(xlim=x.lim.cart)
  } else if (! missing(y.lim.cart)) {
    if (length(y.lim.cart) != 2)
      stop("Invalid format, use y.lim.cart = c(1000,NA) !")
    if (which(is.na(y.lim.cart))==1) {
      y.lim.cart[1] = min(df.melt[,y.id])
    } else if (which(is.na(y.lim.cart))==2) {
      y.lim.cart[2] = max(df.melt[,y.id])
    }
    p <- p + coord_cartesian(ylim=y.lim.cart)
  } 
  
  if (! missing(palette)) {
    if (length(palette) == 1)
      p <- p + scale_colour_brewer(palette=palette)
    else if (length(palette) <= 3)
      p <- p + scale_colour_gradientn(colours=palette)
    else 
      p <- p + scale_fill_manual(values=palette)
  }
  
  # cannot use both legend.col and legend.nrow
  if (legend.col > 1 && legend.nrow == 0)
    p <- p + guides(fill=guide_legend(ncol=legend.col))
  if (legend.nrow > 0 && legend.col == 1)
    p <- p + guides(fill=guide_legend(nrow=legend.nrow,byrow=TRUE))
  
  if (x.lab=="x.id") 
    x.lab = x.id
  if (y.lab=="y.id") 
    y.lab = y.id
  if (!missing(legend.title))
    p <- p + labs(fill=legend.title)
  
  p <- p + theme_bw() + xlab(x.lab) + ylab(y.lab) + ggtitle(title) +
    theme(plot.title = element_text(size = title.size),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank())
  
  if (rotate.x.text) 
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
  
  return(p)
}


#' @details 
#' \code{ggBoxWhiskersPlot} creates box Whiskers plot. 
#' 
#' @param outlier.colour The colour of outliers in box whiskers plot 
#' used for \code{outlier.colour} in \code{\link{geom_boxplot}} in \code{\link{ggplot}}. 
#' Default to alpha("black", 0.3).
#' @param y.lower,y.upper The \code{\link{limits}} of y-axis, which can cut off the data points.
#' @keywords graph
#' @export
#' @examples
#' box.plot <- ggBoxWhiskersPlot(df.melt, x.id="test", y.id="performance")
#' 
#' @rdname ggPlot
ggBoxWhiskersPlot <- function(df.melt, x.id, y.id, fill.id, x.facet.id, y.facet.id, 
                              outlier.colour = alpha("black", 0.3), y.scale="identity", 
                              y.lower=NA, y.upper=NA, x.lim.cart, y.lim.cart,
                              palette, rotate.x.text=FALSE, dodge.width=0.8,
                              title="Box Whiskers Plot", title.size = 10, x.lab="x.id", y.lab="y.id") {
  if (!is.element(tolower(x.id), tolower(colnames(df.melt))))
    stop(paste0("Data frame do NOT have column name \"", x.id, "\" !"))
  if (!is.element(tolower(y.id), tolower(colnames(df.melt))))
    stop(paste0("Data frame do NOT have column name \"", y.id, "\" !"))
  
  require(ggplot2)
  if (missing(fill.id)) {
    p <- ggplot(df.melt, aes_string(x = x.id, y = y.id)) +
      geom_boxplot(outlier.colour = alpha("black", 0.3))
  } else {
    if (!is.element(tolower(fill.id), tolower(colnames(df.melt))))
      stop(paste0("Data frame do NOT have column name \"", fill.id, "\" !"))
    
    p <- ggplot(df.melt, aes_string(x = x.id, y = y.id, fill = fill.id)) +
      geom_boxplot(outlier.colour = alpha("black", 0.3), position=position_dodge(width=dodge.width))
  }
  p <- p + scale_shape(solid = FALSE) 
  #geom_jitter(alpha = 0.5) + 
  
  if (! missing(x.facet.id) || ! missing(y.facet.id)) {
    fac1 <- "."
    fac2 <- "."
    if (! missing(x.facet.id))
      fac1 <- x.facet.id
    if (! missing(y.facet.id))
      fac2 <- y.facet.id
    
    p <- p + facet_grid(reformulate(fac2,fac1))
  }
  
  if (y.scale=="identity") {
    p <- p + scale_y_continuous(expand = c(0,0))
  } else if (y.scale=="per") {
    require(scales)
    p <- p + scale_y_continuous(expand = c(0,0), labels = percent_format()) 
  } else {
    y.breaks <- ComMA::get_breaks_positive_values(max(df.melt[,y.id], start=c(0)))
    p <- p + scale_y_continuous(trans=y.scale, expand = c(0,0), labels = ComMA::scientific_10, breaks = y.breaks) 
  }
  
  if (!is.na(y.lower) || !is.na(y.upper)) 
    p <- p + ylim(y.lower, y.upper) 
  
  if (! missing(x.lim.cart)) {
    if (length(x.lim.cart) != 2)
      stop("Invalid format, use x.lim.cart = c(1000,NA) !")
    if (which(is.na(x.lim.cart))==1) {
      x.lim.cart[1] = min(df.melt[,x.id])
    } else if (which(is.na(x.lim.cart))==2) {
      x.lim.cart[2] = max(df.melt[,x.id])
    }
    p <- p + coord_cartesian(xlim=x.lim.cart)
  } else if (! missing(y.lim.cart)) {
    if (length(y.lim.cart) != 2)
      stop("Invalid format, use y.lim.cart = c(1000,NA) !")
    if (which(is.na(y.lim.cart))==1) {
      y.lim.cart[1] = min(df.melt[,y.id])
    } else if (which(is.na(y.lim.cart))==2) {
      y.lim.cart[2] = max(df.melt[,y.id])
    }
    p <- p + coord_cartesian(ylim=y.lim.cart)
  } 
  
  if (! missing(palette)) {
    if (length(palette) == 1)
      p <- p + scale_colour_brewer(palette=palette)
    else if (length(palette) <= 3)
      p <- p + scale_colour_gradientn(colours=palette)
    else 
      p <- p + scale_fill_manual(values=palette)
  }
  
  if (x.lab=="x.id") 
    x.lab = x.id
  if (y.lab=="y.id") 
    y.lab = y.id
  
  p <- p + theme_bw() + xlab(x.lab) + ylab(y.lab) + ggtitle(title) +
    theme(plot.title = element_text(size = title.size),
          strip.background = element_blank(), panel.grid.minor = element_blank())
  
  if (rotate.x.text) 
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
  
  return(p)
}


#' @details 
#' \code{gtScatterPlot} uses one-line function to plot many types of scatter chart.
#' 
#' @param point.size The size of points from \code{\link{point.size}}. Default to 3.
#' @param text.id Label the data points according \code{text.id} column, 
#' such as "Row.names" column after \code{\link{merge}}.
#' @param text.size,hjust,vjust,alpha Refer to aesthetics from \code{\link{geom_text}}.
#' @param colour.id,shaped.id,linked.id The column name in \code{df.melt} to 
#' define how the data points are coloured, shaped, or linked according their values.
#' @param ellipsed.id The column name in \code{df.melt} to define 
#' how to draw ellipse over data points, which is normally same as 
#' \code{colour.id} to show clusters.
#' @param shapes Manually define the shapes of points. Refer to \code{\link{scale_shape_manual}}, 
#' and \url{http://www.cookbook-r.com/Graphs/Shapes_and_line_types/}.
#' @param scale.colours A vector of colours to use for n-colour gradient in 
#' \code{\link{scale_colour_gradientn}}, such as scale.colours=c("red", "green", "blue").
#' Default not to use.
#' @param xintercept,yintercept,linetype Add horizontal or vertical line. 
#' Refer to \code{\link{geom_hline}} or \code{\link{geom_vline}}.
#' @keywords graph
#' @export
#' @examples 
#' 
#' 
#' @rdname ggPlot
gtScatterPlot <- function(df.melt, x.id, y.id, colour.id, shaped.id, linked.id,
                          text.id, ellipsed.id, scale.colours, palette, shapes, 
                          point.size=3, xintercept, yintercept, linetype=2,
                          text.size = 3, hjust=-0.1, vjust = -0.2, text.alpha = 0.5,
                          title.size = 10, title="Scatter Plot", x.lab="x.id", y.lab="y.id") {
  if (!is.element(tolower(x.id), tolower(colnames(df.melt))))
    stop(paste0("Data frame do NOT have column name \"", x.id, "\" !"))
  if (!is.element(tolower(y.id), tolower(colnames(df.melt))))
    stop(paste0("Data frame do NOT have column name \"", y.id, "\" !"))
  
  col.id <- c()
  if (! missing(colour.id)) {
    if (! colour.id %in% colnames(df.melt))
      stop(paste("Invalid colour.id,", colour.id,  "not exsit in column names !"))
    col.id <- c(col.id, colour.id)
  } 
  if (! missing(linked.id)) { 
    if (! linked.id %in% colnames(df.melt) )
      stop(paste("Invalid linked.id,", linked.id,  "not exsit in column names !"))
    col.id <- c(col.id, linked.id)
  } 
  if (! missing(shaped.id)) { 
    if (! shaped.id %in% colnames(df.melt) )
      stop(paste("Invalid shaped.id,", shaped.id,  "not exsit in column names !"))
    col.id <- c(col.id, shaped.id)
  } 
  
  require(ggplot2)
  if (missing(colour.id)) {
    p <- ggplot(df.melt, aes_string(x = x.id, y = y.id)) 
  } else {
    if (!is.element(tolower(colour.id), tolower(colnames(df.melt))))
      stop(paste0("Data frame do NOT have column name \"", colour.id, "\" !"))
    
    p <- ggplot(df.melt, aes_string(x = x.id, y = y.id, colour=colour.id)) 
    if (! missing(scale.colours))
      p <- p + scale_colour_gradientn(colours = scale.colours) 
  }
  
  if (! missing(shaped.id)) {
    if (!is.element(tolower(shaped.id), tolower(colnames(df.melt))))
      stop(paste0("Data frame do NOT have column name \"", shaped.id, "\" !"))
    
    p <- p + geom_point(aes_string(shape=shaped.id), size = point.size) 
    if (! missing(shapes)) {
      p <- p + scale_shape_manual(values=shapes) 
    } else {
      # The shape palette can deal with a maximum of 6 discrete values
      n_shape <- length(unique(df.melt[,shaped.id]))
      myshape <- seq(1, (1 + n_shape-1))
      p <- p + scale_shape_manual(values=myshape) 
    }
  } else {
    p <- p + geom_point(size = point.size)
  }
  
  if (! missing(linked.id)) {
    # Convex hull http://stackoverflow.com/questions/16428962/convex-hull-ggplot-using-data-tables-in-r
    df.melt <- data.table(df.melt, key = linked.id)
    chull.cmd <- parse(text = paste0('df.melt[, .SD[chull(', x.id, ',', y.id, ')], by = "', linked.id, '"')) 
    # hulls <- pts.mds.dt[, .SD[chull(MDS1, MDS2)], by = linked.id]
    hulls <- eval(chull.cmd)
    
    p <- p + geom_polygon(data = hulls, aes_string(mapping=linked.id), fill = NA, alpha = 0.5)
  }
  
  if (! missing(ellipsed.id)) { # normally ellipsed.id == colour.id to show clusters
    if (! ellipsed.id %in% colnames(df.melt) )
      stop(paste("Invalid ellipsed.id,", ellipsed.id,  "not exsit in column names !"))
    p <- p + stat_ellipse(aes_string(mapping = ellipsed.id), type = "t", linetype = 2)
  }
  
  if (! missing(text.id)) {
    if (!is.element(tolower(text.id), tolower(colnames(df.melt))))
      stop(paste0("Data frame do NOT have column name \"", text.id, "\" !"))
    
    if (! missing(colour.id)) {
      p <- p + geom_text(aes_string(label=text.id, colour=colour.id), 
                         size=text.size, hjust=hjust, vjust=vjust, alpha=text.alpha)
    } else {
      p <- p + geom_text(aes_string(label=text.id), 
                         size=text.size, hjust=hjust, vjust=vjust, alpha=text.alpha)
    }
  }
  
  if (! missing(palette)) {
    if (length(palette) == 1)
      p <- p + scale_colour_brewer(palette=palette)
    else if (length(palette) <= 3)
      p <- p + scale_colour_gradientn(colours=palette)
    else 
      p <- p + scale_fill_manual(values=palette)
  }
  
  if (! missing(xintercept))
    p <- p + geom_vline(xintercept=xintercept,linetype=linetype)
  if (! missing(yintercept))
    p <- p + geom_hline(yintercept=yintercept,linetype=linetype) 
  
  if (x.lab=="x.id") 
    x.lab = x.id
  if (y.lab=="y.id") 
    y.lab = y.id
  
  p <- p + theme_bw() + xlab(x.lab) + ylab(y.lab) + ggtitle(title) +
    theme(plot.title = element_text(size = title.size), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank()) 
  #    guides(shape=guide_legend(nrow=2)) +
  
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
gtLine <- function(df.melt, x.id, y.id, group.id, colour.id, shaped.id, text.id,
                   x.scale="identity", y.scale="identity", shapes, palette, 
                   line.size=0.5, line.alpha=0.75, legend.title, 
                   text.size = 3, hjust=-0.1, vjust = -0.2, text.alpha = 0.5,
                   title="Line", title.size = 10, x.lab="x.id", y.lab="y.id") {
  if (!is.element(tolower(x.id), tolower(colnames(df.melt))))
    stop(paste0("Data frame do NOT have column name \"", x.id, "\" !"))
  if (!is.element(tolower(y.id), tolower(colnames(df.melt))))
    stop(paste0("Data frame do NOT have column name \"", y.id, "\" !"))

  require(ggplot2)
  aes.string <- paste("x =", x.id, "y =", y.id)
  if (! missing(group.id)) {
    if (!is.element(tolower(group.id), tolower(colnames(df.melt))))
      stop(paste0("Data frame do NOT have column name \"", group.id, "\" !"))
    
    aes.string <- paste(aes.string, "group =", group.id)
  } 
  if (! missing(colour.id))  {
    if (!is.element(tolower(colour.id), tolower(colnames(df.melt))))
      stop(paste0("Data frame do NOT have column name \"", colour.id, "\" !"))
    
    aes.string <- paste(aes.string, "colour =", colour.id)
  } 
  p <- p + ggplot(df.melt, aes_string(aes.string)) + geom_line(size=line.size, alpha=0.75) 
  
  if (! missing(shaped.id)) {
    if (!is.element(tolower(shaped.id), tolower(colnames(df.melt))))
      stop(paste0("Data frame do NOT have column name \"", shaped.id, "\" !"))
    
#    if (! missing(colour.id)) 
#      p <- p + geom_point(aes_string(shape=shaped.id, colour=colour.id), size = point.size)
#    else 
      p <- p + geom_point(aes_string(shape=shaped.id), size = point.size)
      
      if (! missing(shapes)) {
        p <- p + scale_shape_manual(values=shapes) 
      } else {
        # The shape palette can deal with a maximum of 6 discrete values
        n_shape <- length(unique(df.melt[,shaped.id]))
        myshape <- seq(1, (1 + n_shape-1))
        p <- p + scale_shape_manual(values=myshape) 
      }
  } 
  
  if (x.scale!="identity") {
    x.breaks <- ComMA::get_breaks_positive_values(max(df.melt[,x.id], start=c(0)))
    p <- p + scale_x_continuous(trans=x.scale, expand = c(0,0), labels = ComMA::scientific_10, breaks = x.breaks) 
  }
  
  if (y.scale!="identity") {
    y.breaks <- ComMA::get_breaks_positive_values(max(df.melt[,y.id], start=c(0)))
    p <- p + scale_y_continuous(trans=y.scale, expand = c(0,0), labels = ComMA::scientific_10, breaks = y.breaks) 
  }
  
  if (! missing(text.id)) {
    if (!is.element(tolower(text.id), tolower(colnames(df.melt))))
      stop(paste0("Data frame do NOT have column name \"", text.id, "\" !"))
    
    if (! missing(colour.id)) {
      p <- p + geom_text(aes_string(label=text.id, colour=colour.id), 
                         size=text.size, hjust=hjust, vjust=vjust, alpha=text.alpha)
    } else {
      p <- p + geom_text(aes_string(label=text.id), 
                         size=text.size, hjust=hjust, vjust=vjust, alpha=text.alpha)
    }
  }
  
  if (! missing(palette)) {
    if (length(palette) == 1)
      p <- p + scale_colour_brewer(palette=palette)
    else if (length(palette) <= 3)
      p <- p + scale_colour_gradientn(colours=palette)
    else 
      p <- p + scale_fill_manual(values=palette)
  }
  
  if (x.lab=="x.id") 
    x.lab = x.id
  if (y.lab=="y.id") 
    y.lab = y.id
  
  p <- p + theme_bw() + xlab(x.lab) + ylab(y.lab) + ggtitle(title) +
    theme(axis.line = element_line(colour = "black"),
          plot.title = element_text(size = title.size),
          panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),
          panel.border = element_blank(), panel.background = element_blank()) 
  
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  
  return(p)
}

