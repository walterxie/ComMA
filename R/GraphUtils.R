# GraphUtils
# Author: Walter Xie
# Accessed on 11 Mar 2016

#' @name pdf
#' @title Create pdf 
#'
#' @details \code{pdfGgplot} creates pdf for ggplot object.
#' 
#' @param gg.plot A \code{\link{ggplot}} object.
#' @param fig.path fig.path The full path of pdf file.
#' @param width, height Refer to \code{\link{pdf}}. Default to width=6, height=6.
#' @param useDingbats Defaults to TRUE, which produces smaller and better output. 
#' Setting this to FALSE can work around font display problems in broken PDF viewers.
#' Refer to \code{\link{pdf}}.
#' @keywords utils
#' @export
#' @examples 
#' pdfGgplot(gg.plot, fig.path="fig.pdf", width=8, height=8)  
#' 
#' @rdname pdf
pdfGgplot <- function(gg.plot, fig.path, width=6, height=6, useDingbats=TRUE) {
  pdf(fig.path, width=width, height=height, useDingbats=useDingbats)	
  print(gg.plot)
  invisible(dev.off()) 
}

#' @details \code{pdfGtable} creates pdf for gtable object.
#' 
#' @param gtable A \code{\link{gtable}} object.
#' @keywords utils
#' @export
#' @examples 
#' pdfGtable(g.table, fig.path="fig.pdf", width=8, height=8)  
#' 
#' @rdname pdf
pdfGtable <- function(g.table, fig.path, width=8, height=8, useDingbats=TRUE) {
  require(grid)
  pdf(fig.path, width=width, height=height, useDingbats=useDingbats)	
  print(grid.draw(g.table))
  invisible(dev.off()) 
}


#' @details \code{grid_arrange_shared_legend} shares a legend 
#' between multiple plots using \code{\link{grid.arrange}}. 
#' Return a \code{\link{gtable}} object.
#' @source \code{grid_arrange_shared_legend} is derived from 
#' \url{http://rpubs.com/sjackman/grid_arrange_shared_legend}
#' 
#' @param ... The list of \code{\link{ggplot}} objects.
#' @keywords utils
#' @export
#' @examples 
#' library(ggplot2); library(grid); library(gridExtra)
#' dsamp <- diamonds[sample(nrow(diamonds), 1000), ]
#' p1 <- qplot(carat, price, data=dsamp, colour=clarity)
#' p2 <- qplot(cut, price, data=dsamp, colour=clarity)
#' p3 <- qplot(color, price, data=dsamp, colour=clarity)
#' p4 <- qplot(depth, price, data=dsamp, colour=clarity)
#' grid_arrange_shared_legend(p1, p2, p3, p4)
#' 
#' @rdname pdf
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  require(ggplot2)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  
  require(gridExtra)
  require(grid)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

#' Defining the scale change.
#' It is mostly used to force bars to start from a lower value than 0 in \pkg{ggplot2} \code{\link{geom_bar}} in R
#' @source \url{http://stackoverflow.com/questions/22295253/force-bars-to-start-from-a-lower-value-than-0-in-ggplot-geom-bar-in-r}
#' 
#' @param base The base of logarithm to use. Default to exp(1).
#' @param from The value to start from. Default to 0.
#' @return
#' The scale.
#' @keywords utils
#' @export
#' @examples 
#' # starts from 1e-2
#' scale_y_continuous(trans = mylog_trans(base=10, from=-2))
mylog_trans <- function(base=exp(1), from=0) {
  trans <- function(x) log(x, base)-from
  inv <- function(x) base^(x+from)
  require(scales)
  trans_new("mylog", trans, inv, log_breaks(base=base), domain = c(base^from, Inf))
}

#' Scientific notation
#' 
#' @param x The number (not string type).
#' @return
#' The mathematical \code{\link{expression}} in scientific notation of a given number.
#' @keywords utils
#' @export
#' @examples 
#' #expression(10^04)
#' scientific_10(10000)
scientific_10 <- function(x) {
  require(scales)
  text=gsub("1e\\+00", "1", scientific_format()(x))
  text=gsub("1e\\+01", "10", text)
  text=gsub("0e\\+00", "0", text)
  text=gsub("1e\\-01", "0.1", text)
  text=gsub("1e\\+", "10^", text)
  parse(text=text)
}

#' Breaks of multiples of 10 for positive values.
#' 
#' @param max.v The max value to create breaks.
#' @param start The vector of values before 10. Default to c(0.1, 1).
#' @return
#' The vector of multiples of 10 used for breaks. 
#' Mostly used with \code{\link{scientific_10}} together.
#' @keywords utils
#' @export
#' @examples 
#' get_breaks_positive_values(68759)
#' [1] 1e-01 1e+00 1e+01 1e+02 1e+03 1e+04 1e+05
#' 
#' scale_y_continuous(trans = "log", labels = ComMA::scientific_10, 
#' breaks = ComMA::get_breaks_positive_values(max(df, start=c(0))))
get_breaks_positive_values <- function(max.v, start=c(0.1, 1)) {
  breaks=c(start, 10, 100, 1000, 10000, 100000, 1000000)
  if (max.v < 50) {
    breaks=c(start, 10)
  } else if (max.v < 500) {
    breaks=c(start, 10, 100)
  } else if (max.v < 5000) {
    breaks=c(start, 10, 100, 1000)
  } else if (max.v < 50000) {
    breaks=c(start, 10, 100, 1000, 10000)
  } else if (max.v < 500000) {
    breaks=c(start, 10, 100, 1000, 10000, 100000)
  } else if (max.v < 5000000) {
    breaks=c(start, 10, 100, 1000, 10000, 100000, 1000000)
  } 
  return(breaks)
}


#' Customised palette with about 70 colours.
#' 
#' @param add.grey Add grey colours in the begining of palette. Default to 0.
#' @param n The number of colors (>= 1) to be in the palette. If missing, then return the whole palette.
#' @return
#' The customised palette with about 70 colours.
#' @keywords utils
#' @export
#' @examples 
#' myPalette <- getMyPalette() 
getMyPalette<-function(n, add.grey=0){
  greyPalette <- c("#CCCCCC", "#999999")
  myPalette <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#AE22EC", "#FAF629",
                 "#CAB2D6", "#6A3D9A", "#B15928", "#E69F00", "#FF6666", "#56B4E9", "#006600", "#009E73", "#F0E442", "#FF3333",    
                 "#0072B2", "#D55E00", "#CC79A7", "#E6AB02", "#FFFF99", "#33ff99", "#fb9a99", "#556B2F", "#fdbf6f", "#cab2d6", 
                 "#377eb8", "#ff9933", "#4daf4a", "#984ea3", "#cc0033", "#a65628", "#9999FF", "#f781bf","#26ED26", "#1f78b4", 
                 "#b2df8a", "#a6cee3", "#e41a1c", "#ffff33", "#99AC06", "#2F3C56", "#A229FA", "#EC22DE", "#1DD3E8", "#9933FF",
                 "#CDC0B0", "#7FFFD4", "#D2691E", "#458B74", "#6495ED", "#FFF8DC", "#00FFFF", "#EE3B3B", "#CAFF70", "#8B008B",
                 "#00008B", "#A52A2A", "#EEC591", "#98F5FF", "#BF3EFF", "#8B4513", "#EE6A50", "#66CD00", "#8B0000", "#e31a1c")
  
  if (add.grey>0 && add.grey<=length(greyPalette))
    myPalette <- c(greyPalette[1:add.grey], myPalette)
  
  if (missing(n) || n < 1 || n > length(myPalette)) 
    n <- length(myPalette)
    
  return(myPalette[1:n])
}

#' Extract legend in \pkg{ggplot2}
#' @source \url{http://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2}
#' 
#' @param a.gplot The \code{\link{ggplot}} object.
#' @return
#' The legend.
#' @keywords utils
#' @export
#' @examples 
#' library(ggplot2); library(grid)
#' my_hist<-ggplot(diamonds, aes(clarity, fill=cut)) + geom_bar() 
#' legend <- g_legend(my_hist) 
#' grid.draw(legend) 
g_legend<-function(a.gplot){
  require(ggplot2)
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
