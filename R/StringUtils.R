# String Utils
# Author: Walter Xie
# Accessed on 15 Sep 2016

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
#' @param y The string of a replacement, default is the empty string.
#' @keywords utils
#' @export
#' @examples 
#' text = "   foo bar  baz 3 "
#' trimStartEndSpace(text)
#' [1] "foo bar  baz 3"
#'
#' @rdname stringUtils
trimStartEndSpace <- function (x, y="") gsub("^\\s+|\\s+$", y, x)

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
trimSpace <- function (x, y="") gsub("[[:space:]]", y, x)
trimSpace.old <- function (x, y="") gsub("\\s", y, x)

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

#' @details \code{getFirstNInString} returns the first \code{n} elements in a string, 
#' using \code{\link{paste}}. 
#' If \code{n} is greater than length,
#' then return the whole vector or list in a string.
#' 
#' @param vect A Vector or List.
#' @param n The first n elements.
#' @export
#' @examples 
#' getFirstNInString(1:10, 3)
#' getFirstNInString(1:10, 20)
#' 
#' @rdname stringUtils
getFirstNInString <- function(vect, n, collapse=", ", tail="...") {
  if (length(vect) < n)
    return( paste(vect, collapse=collapse) )
  else
    return( paste(paste(vect[1:n], collapse=collapse), tail) )
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
#' and their frequency given a vector \code{x}, if \code{x} is a data frame, 
#' then return unique pairs counts including zero counts. 
#' Refere to \url{http://stackoverflow.com/questions/16905425/find-duplicate-values-in-r}.
#' 
#' @export
#' @examples 
#' freqUniqueValue(c("a", "b", "a"))
#' #    Var1 Freq
#' # 1    a    2
#' # 2    b    1
#' freqUniqueValue(data.frame(x=c(1, 1, 2), y=c(3, 4, 3)), rm.zero=F)
#' #   x y Freq
#' # 1 1 3    1
#' # 2 2 3    1
#' # 3 1 4    1
#' # 4 2 4    0
#' 
#' @rdname stringUtils
freqUniqueValue <- function(x, rm.zero=TRUE) {
  freq.df <- data.frame(table(x))
  if (rm.zero) freq.df <- freq.df[freq.df$Freq>0,]
  freq.df
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

#' @details \code{gsubDF} applies \code{\link{gsub}} to 
#' the entire data frame.
#' 
#' @param df Data frame as the input.
#' @param pattern,replacement,... Arguemnts for \code{\link{gsub}}.
#' @export
#' @examples 
#' gsubDF(".00", "", df)
#' 
#' @rdname stringUtils
gsubDF <- function(pattern, replacement, df, stringsAsFactors=FALSE, check.names=FALSE, ...) {
  df1 <- data.frame(apply(df, 2, function(x) gsub(pattern, replacement, as.character(x), ...)), 
                    stringsAsFactors=stringsAsFactors, check.names=check.names)
  rownames(df1) <- rownames(df)
  return(df1)
}

#' @details \code{pasteDF} \code{\link{paste}} strings to 
#' the entire data frame. 
#' This will change the type to \code{\link{character}}.
#' 
#' @param sep,collapse Arguemnts for \code{\link{paste}}.
#' @export
#' @examples 
#' pasteDF(df, "%")
#' 
#' @rdname stringUtils
pasteDF <- function(df, ..., sep = " ", collapse = NULL, stringsAsFactors=FALSE, check.names=FALSE) {
  df1 <- data.frame(apply(df, 2, function(x) paste(as.character(x), ..., sep=sep, collapse=collapse)), 
                    stringsAsFactors=stringsAsFactors, check.names=check.names)
  rownames(df1) <- rownames(df)
  colnames(df1) <- colnames(df)
  return(df1)
}
