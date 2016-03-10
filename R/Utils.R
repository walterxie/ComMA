
# Author: Walter Xie
# Accessed on 10 Mar 2016


# http://stackoverflow.com/questions/22295253/force-bars-to-start-from-a-lower-value-than-0-in-ggplot-geom-bar-in-r
# defining the scale change
# scale_y_continuous(trans = mylog_trans(base=10, from=-2)) # starts from 1e-2
mylog_trans <- function(base=exp(1), from=0) {
  require(scales)
  trans <- function(x) log(x, base)-from
  inv <- function(x) base^(x+from)
  trans_new("mylog", trans, inv, log_breaks(base=base), domain = c(base^from, Inf))
}

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Fix x-axis number format
scientific_10 <- function(x) {
  text=gsub("1e\\+00", "1", scientific_format()(x))
  text=gsub("1e\\+01", "10", text)
  text=gsub("1e\\+", "10^", text)
  parse(text=text)
}

######## get table in the format of "corr (sign)" #######
# Input: corr.sign.matrix is a matrix having same row and col names, 
# lower triangle is correlations (or equivalent), upper triangle is significance
# Output corr.sign.table is in the format of "corr (sign)"
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

# http://stackoverflow.com/questions/2261079/how-to-trim-leading-and-trailing-whitespace-in-r 
# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#http://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r
#Extracting the last n characters from a string x
substrRight <- function(x, n) {substr(x, nchar(x)-n+1, nchar(x))}

# get plural for report
getPlural <- function (singular) {
  if (tolower(singular)=="species") {
    plural = singular
  } else if (tolower(singular)!="phylumn") {
    plural = "phyla"
  } else {
    plural = paste(singular,"s",sep="")
  }
  return(plural) 
}

# http://stackoverflow.com/questions/29214932/split-a-file-path-into-folder-names-vector
split_path <- function(path) {
  rev(setdiff(strsplit(path,"/|\\\\")[[1]], ""))
}

#' Generate coordinates for 2 clusters. 
#' 
#' Generate coordinates for 2 clusters.
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
