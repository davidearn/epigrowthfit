\name{print.egf}
\alias{print.egf}
\title{Printing Model Objects}
\description{
A method for printing objects of class \code{\link{egf}}.
}
\usage{
\method{print}{egf}(x, width = 0.9 * getOption("width"), indent = 2L, \dots)
}
\arguments{
\item{x}{
  an \code{\link{egf}} object.}
\item{width}{
  an integer width for header text.}
\item{indent}{
  an integer indent for body text.}
\item{\dots}{
  unused optional arguments.}
}
\value{
The argument \code{x}, unchanged but \link{invisible}.
}
\seealso{
The generic function \code{\link{print}}.
}
