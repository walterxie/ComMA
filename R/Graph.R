# Graph
# Author: Walter Xie
# Accessed on 10 Mar 2016

#' @name ggPlot
#' @title One-line command to get \pkg{ggplot} 
#'
#' @description Simplify \pkg{ggplot} codes into functions that can get a chart from one-line command.
#' 
#' @details 
#' \code{ggLine} adds a line to the given \code{\link{ggplot}} object.
#' 
#' @param xintercept,yintercept,intercept,slope,smooth.method Refer to \pkg{ggplot2} 
#' \code{\link{geom_vline}}, \code{\link{geom_hline}}, \code{\link{geom_abline}}, 
#' \code{\link{geom_smooth}}. They cannot be used at the same time.
#' @param linetype \code{\link{linetype}} 0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash.
#' @keywords graph
#' @export
#' @examples 
#' p <- ggLine(p, linetype = 2, yintercept = 1)
#' p <- ggLine(p, smooth.method = "lm")
#' 
#' @rdname ggPlot
ggLine <- function(gg.plot, linetype=1, xintercept, yintercept, intercept, slope, smooth.method) {
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
#' \code{ggScatterPlot} uses one-line function to plot many types of scatter chart.
#' 
#' @param limits A numeric vector of length two providing limits of the scale. 
#' Use NA to refer to the existing minimum or maximum. Refer to \code{\link{scale_y_continuous}}.
#' Set lower bound only to y-axis using limits = c(1000, NA). Defaul to NULL. 
#' @keywords graph
#' @export
#' @examples 
#' 
#' 
#' @rdname ggPlot
ggScatterPlot <- function(df.melt, x.id, y.id, fill.id, point.size=3, palette=NULL, 
                          rotate.x.text=FALSE, legend.col=1, legend.nrow=0, limits=NULL,
                          title="Bar Chart", title.size = 10, x.lab="x.id", y.lab="y.id", 
                          addLabel=FALSE, label.size = 3, hjust=-0.1, vjust = -0.2, alpha = 0.5) {
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
  p <- p + geom_point(size=point.size) + #aes(shape=factor(species)),
  
  if (y.scale=="nor") {
    p <- p + scale_y_continuous(expand = c(0,0), limits=limits)
  } else if (y.scale=="per") {
    require(scales)
    p <- p + scale_y_continuous(expand = c(0,0), limits=limits, labels = percent_format()) 
  } else {
    y.breaks <- ComMA::get_breaks_positive_values(max(df.melt[,y.id], start=c(0)))
    p <- p + scale_y_continuous(trans=y.scale, expand = c(0,0), limits=limits, 
                                labels = ComMA::scientific_10, breaks = y.breaks) 
  }
  
  if (!is.null(palette))
    p <- p + scale_fill_manual(values=palette)
  # cannot use both legend.col and legend.nrow
  if (legend.col > 1 && legend.nrow == 0)
    p <- p + guides(fill=guide_legend(ncol=legend.col))
  if (legend.nrow > 0 && legend.col == 1)
    p <- p + guides(fill=guide_legend(nrow=legend.nrow,byrow=TRUE))
  
  if (x.lab=="x.id") 
    x.lab = x.id
  if (y.lab=="y.id") 
    y.lab = y.id
  
  p <- p + theme_bw() + xlab(x.lab) + ylab(y.lab) + ggtitle(title) +
    theme(plot.title = element_text(size = title.size),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank())
  
  if (rotate.x.text) 
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
  
  return(p)
}

#' @details 
#' One-line bar chart function to plot many types of bar chart, 
#' such as normal bars, log-scaled bars, percentage bars, and also grouping.
#' @param df.melt A data frame already \code{\link{melt}}. 
#' @param x.id,y.id, fill.id The string of column names in the data frame 
#' used for \code{x, y, fill} in \code{\link{aes}} in \code{\link{ggplot}}.
#' @param bar.pos Position adjustment for \code{\link{geom_bar}}, either as a string, 
#' or the result of a call to a position adjustment function. Default to "dodge". 
#' Use "fill" to generate group percentage bars.
#' @param y.scale The string defines the data scale used in y-axis, 
#' which can be "nor" standing for normal, or "per" standing for percentage, 
#' moreover either the name of a transformation object for \code{\link{scale_y_continuous}} 
#' (e.g. \code{trans="log"}), or the object itself. Built-in transformations include 
#' "asn", "atanh", "boxcox", "exp", "identity", "log", "log10", "log1p", "log2", "logit", 
#' "probability", "probit", "reciprocal", "reverse" and "sqrt". Default to "nor". 
#' @param x.lim.cart,y.lim.cart Setting limits on the coordinate system will zoom the plot, 
#' and will not change the underlying data like setting limits on a scale will. 
#' Refer to \code{\link{coord_cartesian}}. Set lower bound only to y-axis using y.lim.cart=c(1000,NA). 
#' Defaul to NULL. 
#' @param title Graph title
#' @param x.lab,y.lab The label of x-axis or y-axis, such as plot names.
#' @param palette The colour palette for bars. 
#' If NULL, then use \code{\link{ggplot}} default colours. Default to NULL. 
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
#' @rdname ggPlot
ggBarChart <- function(df.melt, x.id, y.id, fill.id, bar.pos="dodge", y.scale="nor", palette=NULL, 
                       rotate.x.text=FALSE, legend.col=1, legend.nrow=0, x.lim.cart=NULL, y.lim.cart=NULL,
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
  
  if (y.scale=="nor") {
    p <- p + scale_y_continuous(expand = c(0,0))
  } else if (y.scale=="per") {
    require(scales)
    p <- p + scale_y_continuous(expand = c(0,0), labels = percent_format()) 
  } else {
    y.breaks <- ComMA::get_breaks_positive_values(max(df.melt[,y.id], start=c(0)))
    p <- p + scale_y_continuous(trans=y.scale, expand = c(0,0), 
                                labels = ComMA::scientific_10, breaks = y.breaks) 
  }

  if (!is.null(x.lim.cart)) {
    if (which(is.na(x.lim.cart))==1) {
      x.lim.cart[1] = min(df.mu.aggr[,x.id])
    } else if (which(is.na(x.lim.cart))==2) {
      x.lim.cart[2] = max(df.mu.aggr[,x.id])
    }
    p <- p + coord_cartesian(xlim=x.lim.cart)
  } else if (!is.null(y.lim.cart)) {
    if (which(is.na(y.lim.cart))==1) {
      y.lim.cart[1] = min(df.mu.aggr[,y.id])
    } else if (which(is.na(y.lim.cart))==2) {
      y.lim.cart[2] = max(df.mu.aggr[,y.id])
    }
    p <- p + coord_cartesian(ylim=y.lim.cart)
  } 
  
  if (!is.null(palette))
    p <- p + scale_fill_manual(values=palette)
  # cannot use both legend.col and legend.nrow
  if (legend.col > 1 && legend.nrow == 0)
    p <- p + guides(fill=guide_legend(ncol=legend.col))
  if (legend.nrow > 0 && legend.col == 1)
    p <- p + guides(fill=guide_legend(nrow=legend.nrow,byrow=TRUE))
  
  if (x.lab=="x.id") 
    x.lab = x.id
  if (y.lab=="y.id") 
    y.lab = y.id
  
  p <- p + theme_bw() + xlab(x.lab) + ylab(y.lab) + ggtitle(title) +
    theme(plot.title = element_text(size = title.size),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank())
  
  if (rotate.x.text) 
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
  
  return(p)
}

#' @details 
#' Box Whiskers Plot. 
#' 
#' @param outlier.colour The colour of outliers in box whiskers plot 
#' used for \code{outlier.colour} in \code{\link{geom_boxplot}} in \code{\link{ggplot}}. 
#' Default to alpha("black", 0.3).
#' @keywords graph
#' @export
#' @examples
#' box.plot <- ggBoxWhiskersPlot(df.melt, x.id="test", y.id="performance")
#' 
#' @rdname ggPlot
ggBoxWhiskersPlot <- function(df.melt, x.id, y.id, fill.id, outlier.colour = alpha("black", 0.3), y.scale="nor", 
                              y.lower=NA, y.upper=NA, palette=NULL, rotate.x.text=FALSE, dodge.width=0.8,
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
  
  if (!is.null(palette))
    p <- p + scale_fill_manual(values=palette)
  
  if (y.scale=="nor") {
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

#' Percentage bar chart coloured by groups. 
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
  
  p <- ggBarChart(df.melt, x.id="variable", y.id="value", fill.id=melt.id, bar.pos="fill", 
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
#' @param x.lab,y.lab The label of x-axis or y-axis, such as plot names.
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
ggBarYAcrossX <- function(df.aggre, melt.id="samples", print.xtable=NULL, title="The number of OTUs across the number of samples", 
                       x.lab="Number of samples crossed", y.lab="Number of OTUs", 
                       legend.title="", legend.labels=c("OTUs", "reads"), x.lab.interval=1) {
  if (!is.element(tolower(melt.id), tolower(colnames(df.aggre))))
    stop(paste0("Data frame column names do NOT have \"", melt.id, "\" for melt function !"))
  
  require(reshape2)
  df.melt <- melt(df.aggre, id=melt.id)
  
  x.breaks <- seq(min(df.aggre[,melt.id]), max(df.aggre[,melt.id]), x.lab.interval)
  y.breaks <- ComMA::get_breaks_positive_values(max(df.aggre, start=c(0)))
  
  p <- ggBarChart(df.melt, x.id=melt.id, y.id="value", fill.id="variable", 
                  y.scale="log", title=title, x.lab=x.lab, y.lab=y.lab)
  
  # conside x as discrete values
  p <- p + scale_x_discrete(breaks=x.breaks) +
    scale_fill_discrete(legend.title, labels=legend.labels) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black", size = 0.3),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
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
