\name{egf1}
\alias{egf1}
\title{Fit Nonlinear Models of Epidemic Growth}
\description{
An interface to \code{\link{egf}} providing simpler usage for the case
in which there is exactly one time series and exactly one fitting
window.
}
\usage{
egf1(model, time, x, start = NULL, end = NULL,
     init = NULL, map = NULL, \dots)
}
\arguments{
\item{model}{
  an \R{} object of class \code{\link{egf_model}} specifying a top level
  nonlinear model.}
\item{time}{
  an increasing numeric vector.
  \code{\link{Date}}, \code{\link{POSIXct}}, and \code{\link{POSIXlt}}
  vectors are supported and coerced to numeric with
  \code{\link{julian}(time)}.}
\item{x}{
  a non-negative numeric vector of length equal to \code{length(time)}
  such that \code{x[i]} is the number of cases observed over the
  interval \code{(time[i-1], time[i]]}.  \code{x[1]} is ignored
  internally.  Elements that are not integer-valued are rounded with a
  warning.}
\item{start, end}{
  numeric (or \code{\link{Date}}, \code{\link{POSIXct}},
  \code{\link{POSIXlt}}) vectors of length 1 specifying a fitting
  window.  Time points before \code{start} or after \code{end} are
  discarded.  By default, no time points are discarded.}
\item{init}{
  a numeric vector of length equal to the number of nonlinear model
  parameters specifying values to be used in the first likelihood
  evaluation.  The mapping of elements to parameters is defined by
  \code{\link{egf_top}(model)}.  Use \code{NA} to indicate that elements
  should retain their default value computed from \code{time} and
  \code{x}.}
\item{map}{
  an index vector for \code{init}.  The indexed parameters are fixed at
  their initial values.}
\item{\dots}{
  additional arguments passed to \code{\link{egf}}.}
}
\value{
A list inheriting from class \code{egf}.
See topic \code{\link{egf-class}} for class documentation.
}
\seealso{
The many methods for class \code{egf},
listed by \code{\link{methods}(class = "egf")}.

Function \code{\link{egf}} for handling collections of time series.
}
\examples{
## TODO
}
