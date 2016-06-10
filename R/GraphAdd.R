# Graph
# Author: Walter Xie
# Accessed on 17 Apr 2016

#' Add a line
#' 
#' Add a line \code{\link{geom_line}} to a given \code{\link{ggplot}} object.
#' 
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
#' p <- ggHeatmap(ranks.by.group, melt.id="plot")
#' p <- p + ggAddLine(linetype = 2, yintercept = 1)
#' p <- p + ggAddLine(smooth.method = "lm")
ggAddLine <- function(linetype=1, xintercept, yintercept, intercept, slope, smooth.method, ...) {
  if (!missing(xintercept)) {
    geom_vline(linetype=linetype, xintercept = xintercept, ...)
  } else if (!missing(yintercept)) {
    geom_hline(linetype=linetype, yintercept = yintercept, ...)
  } else if (!missing(intercept) && !missing(slope)){
    geom_abline(linetype=linetype, intercept = intercept, slope = slope, ...)
  } else if (!missing(smooth.method)){  
    geom_smooth(linetype=linetype, method = smooth.method, se = FALSE, ...)
  } else {
    stop("Invalid input !")
  }
}

#' Add numbers
#' 
#' Add numbers as text in a \code{\link{ggplot}} object, such as mean of box plot. 
#' Refer to \code{\link{stat_summary}}.
#' 
#' @param label.id A string of a data mapping (column name) provides labels. 
#' If NULL, the default, then use \code{\link{stat_summary}} by given \code{fun.y.lab}. 
#' Othewise, use \code{\link{geom_text}} by given \code{label.id}.
#' @param hjust,vjust The adjustment of label positions. Only activated if label.id is not NULL.
#' @param dodge.width Dodging width, when different to the width of the individual elements. 
#' Default to 0.8. Refer to \code{\link{position_dodge}}.
#' @param text.size The text size of labels. Default to 3.
#' @param text.colour The text colour. Default to black.
#' @param ... Other arguments passed to \code{\link{geom_text}} or \code{\link{stat_summary}}.
#' @keywords graph
#' @export
#' @examples 
#' p <- ggBarChart(myData, x.id="model", y.id="mean", fill.id="OS", y.lab="performance")
#' p <- p + ggAddNumbers(fun.y.lab=length, y.adj=1.02)
ggAddNumbers <- function(label.id=NULL, hjust=0, vjust=0, 
                         dodge.width=0.9, text.size=4, text.colour="black", ...) {
#    if (! label.id %in% mapping)
#      stop("Cannot find label.id ", label.id, " from data$mapping in ggplot !")
    geom_text(aes_string(label=label.id, hjust=hjust, vjust=vjust), position = position_dodge(width=dodge.width), 
                             colour = text.colour, size = text.size, ...)
}

#' Add numbers
#' 
#' Add numbers as text in a \code{\link{ggplot}} object, such as mean of box plot. 
#' Refer to \code{\link{stat_summary}}.
#' 
#' @param fun.y.lab A function to calculate numbers displayed in the figure. 
#' Only activated if label.id is NULL. Default to function \code{\link{mean}}. 
#' Ues \code{\link{length}} to show number of observations.
#' @param fun.y.pos A function to calculate the initial position of text on y-value. 
#' Only activated if label.id is NULL. Default to \code{\link{median}}.
#' @param y.adj The propotion of the initial position of text on y-value. 
#' Only activated if label.id is NULL.
#' > 1 will raises the text, and < 1 will sinks the text. Default to 0.98.
#' @param digits Integer indicating the number of decimal places for \code{\link{round}}.
#' Only activated if label.id is NULL.
#' @param dodge.width Dodging width, when different to the width of the individual elements. 
#' Default to 0.8. Refer to \code{\link{position_dodge}}.
#' @param text.size The text size of labels. Default to 3.
#' @param text.colour The text colour. Default to black.
#' @param ... Other arguments passed to \code{\link{geom_text}} or \code{\link{stat_summary}}.
#' @keywords graph
#' @export
#' @examples 
#' bar.chart + ggAddNumbers(fun.y.lab=mean, y.adj=0.9)
# # "total" have to exist in p$mapping
# p <- p + ggAddNumbersFun(label.id="total", hjust=ifelse(sign(total)>0, 1, 0))
ggAddNumbersFun <- function(fun.y.lab=mean, fun.y.pos=median, y.adj=0.98, digits=2, 
                            dodge.width=0.9, text.size=4, text.colour="black", ...) {
    stat_summary(fun.data = function(y) {return( c(y = fun.y.pos(y)*y.adj, label = round(fun.y.lab(y),digits)) )}, 
                 geom = "text", position = position_dodge(width=dodge.width), 
                 colour = text.colour, size = text.size, ...)
}


#' Add error bars
#' 
#' Add error bars to a \code{\link{ggplot}} object, mostly points \code{\link{geom_point}}. 
#' Refer to \code{\link{geom_errorbar}}.
#' 
#' @param lower,upper The vector of lower or upper plotted in error bars. May use formula 
#' \code{lower = mean - standard error of mean}, \code{upper = mean + standard error of mean} 
#' to calculate. \code{standard error of mean = sd(x) / sqrt(length(x))}. 
#' @param err.bar.width The width input for \code{\link{geom_errorbar}}.
#' @param dodge.width Dodging width, when different to the width of the individual elements. 
#' Default to 0.9. Refer to \code{\link{position_dodge}}.
#' @param ... Other arguments passed to \code{\link{geom_errorbar}}.
#' @keywords graph
#' @export
#' @examples 
#' # lower = mean - standard error of mean, upper = mean + standard error of mean
#' p <- p + ggAddErrorBars(lower, upper)
#' 
#' # Replicate of http://www.r-bloggers.com/building-barplots-with-error-bars/
#' perf.df <- ComMA::readFile("data-raw/model.test.txt")
#' myData <- aggregate(perf.df$performance, by = list(model = perf.df$model, OS = perf.df$OS),
#'     FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
#' myData <- do.call(data.frame, myData)
#' myData$se <- myData$x.sd / sqrt(myData$x.n)
#' colnames(myData) <- c("model", "OS", "mean", "sd", "n", "se")
#' bar.chart <- ggBarChart(myData, x.id="model", y.id="mean", fill.id="OS", y.lab="performance")
#' lower=myData$mean - myData$se
#' upper=myData$mean + myData$se
#' bar.chart <- bar.chart + ggAddErrorBars(lower=lower, upper=upper)
#' bar.chart
ggAddErrorBars <- function(lower, upper, dodge.width=0.9, err.bar.width=0.25, ...) {
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width=dodge.width), 
                width = err.bar.width, ...)
}