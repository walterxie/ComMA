% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utils.R
\name{firstN}
\alias{firstN}
\title{Return the first \code{n} elements. 
If \code{n} is greater than length,
then return the whole vector or list.}
\usage{
firstN(vect, n)
}
\arguments{
\item{vect}{A Vector or List.}

\item{n}{The first n elements.}
}
\description{
Return the first \code{n} elements. 
If \code{n} is greater than length,
then return the whole vector or list.
}
\examples{
first.n <- firstN(1:10, 3)
firstN(1:10, 20)
}
\keyword{utils}

