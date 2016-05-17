# Graph
# Author: Walter Xie
# Accessed on 17 Apr 2016

#' Add a line
#' 
#' Add a line \code{\link{geom_line}} to a given \code{\link{ggplot}} object.
#' 
#' @param gg.plot A \code{\link{ggplot}} object.
#' @param xintercept,yintercept,intercept,slope,smooth.method 
#' Refer to \code{\link{geom_vline}}, \code{\link{geom_hline}}, 
#' \code{\link{geom_abline}}, \code{\link{geom_smooth}}. 
#' They cannot be used at the same time.
#' @param linetype Refer to \code{\link{linetype}}: 
#' 0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 
#' 4 = dotdash, 5 = longdash, 6 = twodash.
#' @param ... Other arguments passed to \code{\link{layer}}.
#' @keywords graph
#' @export
#' @examples 
#' p <- ggAddLine(p, linetype = 2, yintercept = 1)
#' p <- ggAddLine(p, smooth.method = "lm")
ggAddLine <- function(gg.plot, linetype=1, xintercept, yintercept, intercept, slope, smooth.method, ...) {
  if (!missing(xintercept)) {
    p <- gg.plot + geom_vline(linetype=linetype, xintercept = xintercept, ...)
  } else if (!missing(yintercept)) {
    p <- gg.plot + geom_hline(linetype=linetype, yintercept = yintercept, ...)
  } else if (!missing(intercept) && !missing(slope)){
    p <- gg.plot + geom_abline(linetype=linetype, intercept = intercept, slope = slope, ...)
  } else if (!missing(smooth.method)){  
    p <- p + geom_smooth(linetype=linetype, method = smooth.method, se = FALSE, ...)
  } else {
    stop("Invalid input !")
  }
  return(p)
}

#' Add numbers
#' 
#' Add numbers as text in a \code{\link{ggplot}} object, such as mean of box plot. 
#' Refer to \code{\link{stat_summary}}.
#' 
#' @param gg.plot A \code{\link{ggplot}} object.
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
ggAddNumbers <- function(gg.plot, fun.y.lab=mean, fun.y.pos=median, y.adj=0.98, digits=2, 
                         dodge.width=0.8, text.size=3, text.colour="black") {
  p <- gg.plot + stat_summary(fun.data = function(y) {return( c(y = fun.y.pos(y)*y.adj, label = round(fun.y.lab(y),digits)) )}, 
                              geom = "text", position = position_dodge(width=dodge.width), colour = text.colour, size = text.size)
}

