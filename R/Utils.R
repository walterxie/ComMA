# Utils
# Author: Walter Xie
# Accessed on 11 Mar 2016

#' Move rows to the last rows of data frame given a column value 
#' 
#' @param data A data frame.
#' @param col The column, which can be a string (column name) 
#' or number (index of the column).  
#' @param regex Regular expression for \code{\link{grep}} to indetify rows. 
#' Default to "unclassified". Use "^string$" to perform exact match.
#' @param ignore.case Refer to \code{\link{grep}}.
#' @return
#' The data frame of rows moved to the last.
#' @keywords utils
#' @export
#' @examples 
#' 
#' # move all rows having "unclassified" phyla to the last rows
#' data <- mvRowsToLast(data, "phylum")
#' # move all rows having exctly mathced "Bacteria" in column 1 to the last rows
#' data <- mvRowsToLast(data, 1, regex="^Bacteria$")
mvRowsToLast <- function(data, col, regex="unclassified", ignore.case = TRUE) {
  id.match <- grep(regex, data[,col], ignore.case = ignore.case)
  if (length(id.match) > 0)
    data <- data[c(setdiff(1:nrow(data), id.match),id.match),]
  return(data)
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
#' \code{trimStartEndSpace} trims the leading and trailing whitespace in a string.
#' Refer to \url{http://stackoverflow.com/questions/2261079/how-to-trim-leading-and-trailing-whitespace-in-r}.
#'  
#' @param x The string.
#' @keywords utils
#' @export
#' @examples 
#' text = "   foo bar  baz 3 "
#' trimStartEndSpace(text)
#' [1] "foo bar  baz 3"
#'
#' @rdname stringUtils
trimStartEndSpace <- function (x) gsub("^\\s+|\\s+$", "", x)

#' @details \code{trimSpace} removes all whitespace from a string.
#' Refer to \url{http://stackoverflow.com/questions/5992082/how-to-remove-all-whitespace-from-a-string}.
#' 
#' @keywords utils
#' @export
#' @examples 
#' text = "   foo bar  baz 3 "
#' trimSpace(text)
#' [1] "foobarbaz3"
#'
#' @rdname stringUtils
trimSpace <- function (x) gsub("[[:space:]]", "", x)
trimSpace.old <- function (x) gsub("\\s", "", x)

#' @details \code{substrLast} extracts the last n characters from a string x.
#' Refer to \url{http://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r}.
#' 
#' @param n The last n characters
#' @keywords utils
#' @export
#' @examples 
#' x <- "some text in a string"
#' substrLast(x, 6)
#' [1] "string"
#' substrLast(x, 8)
#' [1] "a string"
#'
#' @rdname stringUtils
substrLast <- function(x, n) {substr(x, nchar(x)-n+1, nchar(x))}

#' @details \code{substrLastSplit} extracts the last characters from a string 
#' \code{x} splitted by \code{split}. It can be used to get the file name 
#' from a file path, as default \code{split='/'}. 
#' 
#' @param split character vector (or object which can be coerced to such) containing regular expression(s) 
#' (unless fixed = TRUE) to use for splitting. Detail to \code{\link{strsplit}}.
#' @param ... More parameters passed to \code{\link{strsplit}} from \code{substrLastSplit}.
#' @keywords utils
#' @export
#' @examples 
#' filenames <- list.files(".", pattern="*.java", full.names=TRUE, recursive = TRUE)
#' substrLastSplit(filenames)
#' [1] "BEASTInterfaceTest.java"       "BooleanParameterListTest.java" "IntegerParameterListTest.java"
#'
#' @rdname stringUtils
substrLastSplit <- function(x, split="/", ...) {
  str.splits <- strsplit(x, split=split, ...)
  last.split <- c()
  for (str.spl in str.splits) 
    last.split <- c(last.split, str.spl[[length(str.spl)]])
  return(last.split)
}

#' @details \code{substrBetween} extracts the substring between two given regular expressions.
#' Refer to \url{http://stackoverflow.com/questions/14146362/regex-extract-string-between}.
#' 
#' @param l The regex on the left.
#' @param r The regex on the right.
#' @export
#' @examples 
#' x <- "command took 0:1:34.67 (94.67s total)"
#' substrBetween(x, "\\(", "s total\\)")
#' [1] "string"
#'
#' @rdname stringUtils
substrBetween <- function(x, l, r) {gsub(paste0("^.*", l, "(.*)", r, ".*$"), "\\1", x)}

#' @details \code{simpleCap} capitalizes the first letter of a word string.
#' Refer to \url{http://stackoverflow.com/questions/6364783/capitalize-the-first-letter-of-both-words-in-a-two-word-string}.
#' 
#' @export
#' @examples 
#' simpleCap("BACTERIA")
#' # for multi-words
#' taxaGroups <- c("BACTERIA", "FUNGI", "PROTISTS", "ANIMALIA")
#' sapply(taxaGroups, simpleCap)
#' #BACTERIA      FUNGI   PROTISTS   ANIMALIA 
#' #"Bacteria"    "Fungi" "Protists" "Animalia"
#'
#' @rdname stringUtils
simpleCap <- function(x) {
  s <- strsplit(tolower(x), " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
}

#' @details \code{countCharOccurrences} returns the occurrences of given charater 
#' \code{char} in string \code{s}. Refere to 
#' \url{https://techoverflow.net/blog/2012/11/10/r-count-occurrences-of-character-in-string/}.
#' 
#' @export
#' @examples 
#' countCharOccurrences("a", "application")
#' 
#' @rdname stringUtils
countCharOccurrences <- function(char, s) {
  s2 <- gsub(char,"",s)
  return (nchar(s) - nchar(s2))
}

#' @details \code{freqUniqueValue} returns 2-column data frame of the unique values 
#' and their frequency given a vector \code{x}. Refere to 
#' \url{http://stackoverflow.com/questions/16905425/find-duplicate-values-in-r}.
#' 
#' @export
#' @examples 
#' freqUniqueValue(c("a", "b", "a"))
#' #    Var1 Freq
#' # 1    a    2
#' # 2    b    1
#' freqUniqueValue(c(1, 1, 2, 1))
#' 
#' @rdname stringUtils
freqUniqueValue <- function(x) {
  data.frame(table(x))
}

#' @details \code{findDuplicates} uses \code{freqUniqueValue} to 
#' find the duplicates in a vector \code{x}.
#' 
#' @export
#' @examples 
#' findDuplicates(c("a", "b", "a"))
#' 
#' @rdname stringUtils
findDuplicates <- function(x) {
  n_occur <- freqUniqueValue(x)
  n_occur[n_occur$Freq > 1,1]
}

#' @details \code{gusbDF} applies \code{\link{gsub}} to 
#' the entire data frame.
#' 
#' @export
#' @examples 
#' gusbDF(".00", "", df)
#' 
#' @rdname stringUtils
gusbDF <- function(pattern, replacement, df, stringsAsFactors=FALSE, check.names=FALSE, ...) {
  df1 <- data.frame(lapply(df, function(x) gsub(pattern, replacement, x, ...)), 
                    stringsAsFactors=stringsAsFactors, check.names=check.names)
  rownames(df1) <- rownames(df)
  return(df1)
}

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
getPlural <- function (...) {
  plurals <- c()
  for (singular in list(...)) {
    plurals <- c(plurals, 
                switch(tolower(singular),
                       species = singular,
                       phylumn = "phyla",
                       paste0(singular,"s") ))
  }
  return(plurals) 
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

#' normalize given vector. 
#' 
#' @source \url{https://stat.ethz.ch/pipermail/r-help//2012-October/336676.html}.
#' 
#' @param vec A vector.
#' @keywords utils
#' @export
#' @examples 
#' normed <- as.data.frame(lapply(mtcars, normalize))
#' lapply(normed, range)
normalize <- function(vec) {
  (vec - min(vec, na.rm=TRUE))/(max(vec,na.rm=TRUE) - min(vec, na.rm=TRUE))
}


#' Convert data frame columns to different type, 
#' default to numeric type for columns. 
#' 
#' @source Modified from 
#' \url{http://stackoverflow.com/questions/2288485/how-to-convert-a-data-frame-column-to-numeric-type}.
#' 
#' @param df A data frame.
#' @param MARGIN,FUN Refer to \code{\link{apply}}.
#' @keywords utils
#' @export
#' @examples 
#' df <- convertType(df)
#' df <- convertType(df, FUN=as.character)
convertType <- function(df, FUN=as.numeric, stringsAsFactors=FALSE, check.names=FALSE) {
  df1 <- data.frame(lapply(df, FUN), stringsAsFactors=stringsAsFactors, check.names=check.names)
  rownames(df1) <- rownames(df)
  return(df1)
}
