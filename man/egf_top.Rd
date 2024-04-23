\name{egf_top}
\alias{egf_top}
\alias{egf_top.egf}
\alias{egf_top.egf_no_fit}
\alias{egf_top.egf_model}
\title{Top Level Nonlinear Model Parameter Names}
\description{
Retrieves the names used internally for top level nonlinear model
parameters.
}
\usage{
egf_top(object, \dots)

\method{egf_top}{egf_model}(object, link = TRUE, \dots)
\method{egf_top}{egf}(object, link = TRUE, \dots)
}
\arguments{
\item{object}{
  an \R{} object specifying a top level nonlinear model,
  typically of class \code{\link{egf_model}} or \code{\link{egf}}.}
\item{link}{
  a logical flag.  If \code{TRUE}, then \code{"<link>(<name>)"}
  is returned instead of \code{"<name>"}.}
\item{\dots}{
  unused optional arguments.}
}
\value{
A character vector listing names relevant to \code{object}.
}
