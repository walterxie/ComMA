% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StringUtils.R
\name{stringUtils}
\alias{countCharOccurrences}
\alias{findDuplicates}
\alias{freqUniqueValue}
\alias{getFirstNInString}
\alias{gsubDF}
\alias{pasteDF}
\alias{simpleCap}
\alias{stringUtils}
\alias{substrBetween}
\alias{substrLast}
\alias{substrLastSplit}
\alias{trimSpace}
\alias{trimStartEndSpace}
\title{Utils to deal with strings and characters}
\usage{
trimStartEndSpace(x, y = "")

trimSpace(x, y = "")

substrLast(x, n)

substrLastSplit(x, split = "/", ...)

substrBetween(x, l, r)

simpleCap(x)

getFirstNInString(vect, n, collapse = ", ", tail = "...")

countCharOccurrences(char, s)

freqUniqueValue(x, rm.zero = TRUE)

findDuplicates(x)

gsubDF(pattern, replacement, df, stringsAsFactors = FALSE,
  check.names = FALSE, ...)

pasteDF(df, ..., sep = " ", collapse = NULL, stringsAsFactors = FALSE,
  check.names = FALSE)
}
\arguments{
\item{x}{The string.}

\item{y}{The string of a replacement, default is the empty string.}

\item{n}{The last n characters}

\item{split}{character vector (or object which can be coerced to such) containing regular expression(s) 
(unless fixed = TRUE) to use for splitting. Detail to \code{\link{strsplit}}.}

\item{...}{More parameters passed to \code{\link{strsplit}} from \code{substrLastSplit}.}

\item{l}{The regex on the left.}

\item{r}{The regex on the right.}

\item{vect}{A Vector or List.}

\item{pattern, replacement, ...}{Arguemnts for \code{\link{gsub}}.}

\item{df}{Data frame as the input.}

\item{sep, collapse}{Arguemnts for \code{\link{paste}}.}

\item{n}{The first n elements.}
}
\description{
Useful function to manipulate strings and characters in \code{\link{ComMA}}.
}
\details{
\code{trimStartEndSpace} trims the leading and trailing whitespace in a string.
Refer to \url{http://stackoverflow.com/questions/2261079/how-to-trim-leading-and-trailing-whitespace-in-r}.

\code{trimSpace} removes all whitespace from a string.
Refer to \url{http://stackoverflow.com/questions/5992082/how-to-remove-all-whitespace-from-a-string}.

\code{substrLast} extracts the last n characters from a string x.
Refer to \url{http://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r}.

\code{substrLastSplit} extracts the last characters from a string 
\code{x} splitted by \code{split}. It can be used to get the file name 
from a file path, as default \code{split='/'}.

\code{substrBetween} extracts the substring between two given regular expressions.
Refer to \url{http://stackoverflow.com/questions/14146362/regex-extract-string-between}.

\code{simpleCap} capitalizes the first letter of a word string.
Refer to \url{http://stackoverflow.com/questions/6364783/capitalize-the-first-letter-of-both-words-in-a-two-word-string}.

\code{getFirstNInString} returns the first \code{n} elements in a string, 
using \code{\link{paste}}. 
If \code{n} is greater than length,
then return the whole vector or list in a string.

\code{countCharOccurrences} returns the occurrences of given charater 
\code{char} in string \code{s}. Refere to 
\url{https://techoverflow.net/blog/2012/11/10/r-count-occurrences-of-character-in-string/}.

\code{freqUniqueValue} returns 2-column data frame of the unique values 
and their frequency given a vector \code{x}, if \code{x} is a data frame, 
then return unique pairs counts including zero counts. 
Refere to \url{http://stackoverflow.com/questions/16905425/find-duplicate-values-in-r}.

\code{findDuplicates} uses \code{freqUniqueValue} to 
find the duplicates in a vector \code{x}.

\code{gsubDF} applies \code{\link{gsub}} to 
the entire data frame.

\code{pasteDF} \code{\link{paste}} strings to 
the entire data frame. 
This will change the type to \code{\link{character}}.
}
\examples{
text = "   foo bar  baz 3 "
trimStartEndSpace(text)
[1] "foo bar  baz 3"

text = "   foo bar  baz 3 "
trimSpace(text)
[1] "foobarbaz3"

x <- "some text in a string"
substrLast(x, 6)
[1] "string"
substrLast(x, 8)
[1] "a string"

filenames <- list.files(".", pattern="*.java", full.names=TRUE, recursive = TRUE)
substrLastSplit(filenames)
[1] "BEASTInterfaceTest.java"       "BooleanParameterListTest.java" "IntegerParameterListTest.java"

x <- "command took 0:1:34.67 (94.67s total)"
substrBetween(x, "\\\\(", "s total\\\\)")
[1] "string"

simpleCap("BACTERIA")
# for multi-words
taxaGroups <- c("BACTERIA", "FUNGI", "PROTISTS", "ANIMALIA")
sapply(taxaGroups, simpleCap)
#BACTERIA      FUNGI   PROTISTS   ANIMALIA 
#"Bacteria"    "Fungi" "Protists" "Animalia"

getFirstNInString(1:10, 3)
getFirstNInString(1:10, 20)

countCharOccurrences("a", "application")

freqUniqueValue(c("a", "b", "a"))
#    Var1 Freq
# 1    a    2
# 2    b    1
freqUniqueValue(data.frame(x=c(1, 1, 2), y=c(3, 4, 3)), rm.zero=F)
#   x y Freq
# 1 1 3    1
# 2 2 3    1
# 3 1 4    1
# 4 2 4    0

findDuplicates(c("a", "b", "a"))

gsubDF(".00", "", df)

pasteDF(df, "\%")

}
\keyword{utils}

