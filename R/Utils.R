# Utils
# Author: Walter Xie
# Accessed on 11 Mar 2016


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

#' Get a table in the format of 'corr (sign)'
#' 
#' @param corr.sign.matrix The \code{\link{matrix}} of pairwised correlations and significance, 
#' where the upper triangle is significance and the lower triangle is correlation (or equivalent).
#' @param digits The number of digits to keep.  
#' @return
#' The matrix of strings in the format of 'corr (sign)' for report.
#' @keywords utils
#' @export
#' @examples 
#' 
#' getCorrSignTable(corr.sign.matrix, digits=2)
getCorrSignTable <- function(corr.sign.matrix, digits=3) {
	m.corr <- corr.sign.matrix
	m.corr[upper.tri(m.corr)] <- 0
	m.corr <- formatC(signif(m.corr,digits=digits), digits=digits,format="fg", flag="#")
	m.sign <- t(corr.sign.matrix)
	m.sign[upper.tri(m.sign)] <- 0
	m.sign <- formatC(signif(m.sign,digits=digits), digits=digits,format="fg", flag="#")

	corr.sign.table <- matrix( paste(m.corr, " (", m.sign, ")", sep=""), nrow=nrow(m.corr), dimnames=dimnames(m.corr) )
	corr.sign.table[corr.sign.table=="0 (0)"] <- ""

	corr.sign.table <- corr.sign.table[-1,-ncol(corr.sign.table)]
}

#' @name stringUtils
#' @title Utils to deal with strings and characters
#'
#' @description Useful function to manipulate strings and characters in \code{\link{ComMA}}.
#' 
#' @details 
#' \code{trimStartEnd} trims the leading and trailing whitespace in a string.
#' 
#' @param x The string.
#' @source \url{http://stackoverflow.com/questions/2261079/how-to-trim-leading-and-trailing-whitespace-in-r} 
#' @return
#' The processed string.
#' @keywords utils
#' @export
#' @examples 
#' text = "   foo bar  baz 3 "
#' trimStartEnd(text)
#' [1] "foo bar  baz 3"
#'
#' @rdname stringUtils
trimStartEnd <- function (x) gsub("^\\s+|\\s+$", "", x)

#' @details \code{trimAll} removes all whitespace from a string.
#' 
#' @source \url{http://stackoverflow.com/questions/5992082/how-to-remove-all-whitespace-from-a-string}
#' @keywords utils
#' @export
#' @examples 
#' text = "   foo bar  baz 3 "
#' trimAll(text)
#' [1] "foobarbaz3"
#'
#' @rdname stringUtils
trimAll <- function (x) gsub("\\s", "", x)

#' @details \code{substrRight} extracts the last n characters from a string x
#' 
#' @param n The last n characters
#' @source \url{http://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r}
#' @keywords utils
#' @export
#' @examples 
#' x <- "some text in a string"
#' substrRight(x, 6)
#' [1] "string"
#' substrRight(x, 8)
#' [1] "a string"
#'
#' @rdname stringUtils
substrRight <- function(x, n) {substr(x, nchar(x)-n+1, nchar(x))}

#' @details \code{substrBetween} extracts the substring between two given regular expressions.
#' 
#' @param l The regex on the left.
#' @param r The regex on the right.
#' @source \url{http://stackoverflow.com/questions/14146362/regex-extract-string-between}
#' @examples 
#' x <- "command took 0:1:34.67 (94.67s total)"
#' substrBetween(x, "\\(", "s total\\)")
#' [1] "string"
#'
#' @rdname stringUtils
substrBetween <- function(x, l, r) {gsub(paste0("^.*", l, "(.*)", r, ".*$"), "\\1", x)}

#' Get plural by adding 's' except special plural.
#' 
#' @param singular The singular of a word in string.
#' @return 
#' The plural of given singular. 
#' But it needs to add special plural manully in the code.
#' @keywords utils
#' @export
#' @examples 
#' getPlural("species")
#' [1] "species"
#' getPlural("phylumn")
#' [1] "phyla"
getPlural <- function (singular) {
  if (tolower(singular)=="species") {
    plural = singular
  } else if (tolower(singular)=="phylumn") {
    plural = "phyla"
  } else {
    plural = paste(singular,"s",sep="")
  }
  return(plural) 
}

#' Generate coordinates for 2 clusters. 
#' 
#' @source \url{http://stackoverflow.com/questions/2397097/how-can-a-data-ellipse-be-superimposed-on-a-ggplot2-scatterplot}.
#' 
#' @param n The number of points. Default to 100.
#' @param seed An integer seed for \code{\link{set.seed}}. Default to 101.
#' @keywords utils
#' @export
#' @examples 
#' df.clusters <- random2Clusters()
random2Clusters <- function(n=100, seed=101) {
  #bootstrap
  set.seed(seed)
  x <- rnorm(n, mean=2)
  y <- 1.5 + 0.4*x + rnorm(n)
  df <- data.frame(x=x, y=y, group="A")
  x <- rnorm(n, mean=2)
  y <- 1.5*x + 0.4 + rnorm(n)
  df <- rbind(df, data.frame(x=x, y=y, group="B"))
}

