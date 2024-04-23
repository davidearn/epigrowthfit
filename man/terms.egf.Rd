\name{terms.egf}
\alias{terms.egf}
\alias{terms.egf_no_fit}
\title{Model Terms}
\description{
Extracts the \code{\link[=terms.object]{terms}} object
corresponding to a top level nonlinear model parameter.
}
\usage{
\method{terms}{egf}(x, top = egf_top(x), \dots)
}
\arguments{
\item{x}{
  an \code{\link{egf}} object.}
\item{top}{
  a character string specifying a top level nonlinear model parameter.}
\item{\dots}{
  unused optional arguments.}
}
\value{
A \code{\link[=terms.object]{terms}} object.
}
\seealso{
The generic function \code{\link{terms}}.
}
