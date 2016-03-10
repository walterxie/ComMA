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
  require(scales)
  trans <- function(x) log(x, base)-from
  inv <- function(x) base^(x+from)
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
#' scientific_10(10000)
#' expression(10^04)
scientific_10 <- function(x) {
  require(scales)
  text=gsub("1e\\+00", "1", scientific_format()(x))
  text=gsub("1e\\+01", "10", text)
  text=gsub("1e\\+", "10^", text)
  parse(text=text)
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

#' Trim leading and trailing whitespace in a string.
#' 
#' @param x The string.
#' @source \url{http://stackoverflow.com/questions/2261079/how-to-trim-leading-and-trailing-whitespace-in-r} 
#' @return
#' The string to remove leading or trailing whitespace.
#' @keywords utils
#' @export
#' @examples 
#' text = "   foo bar  baz 3 "
#' trimStartEnd(text)
#' [1] "foo bar  baz 3"
trimStartEnd <- function (x) gsub("^\\s+|\\s+$", "", x)

#' Extracting the last n characters from a string x
#' 
#' @param x The string.
#' @param n The last n characters
#' @source \url{http://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r}
#' @return
#' 
#' @keywords utils
#' @export
#' @examples 
#' x <- "some text in a string"
#' substrRight(x, 6)
#' [1] "string"
#' substrRight(x, 8)
#' [1] "a string"
substrRight <- function(x, n) {substr(x, nchar(x)-n+1, nchar(x))}

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
