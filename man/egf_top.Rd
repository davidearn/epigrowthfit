\name{egf_top}
\alias{egf_top}
\alias{egf_top.default}
\alias{egf_top.egf_model}
\alias{egf_top.egf}
\alias{egf_top.egf_no_fit}
\title{Get top level nonlinear model parameter names}
\usage{
egf_top(object, ...)

\method{egf_top}{default}(object, link = TRUE, ...)

\method{egf_top}{egf_model}(object, link = TRUE, ...)

\method{egf_top}{egf}(object, link = TRUE, ...)

\method{egf_top}{egf_no_fit}(object, link = TRUE, ...)
}
\arguments{
\item{object}{An \R object specifying a top level nonlinear model, or \code{NULL}.}

\item{...}{Unused optional arguments.}

\item{link}{A \link{logical} flag. If \code{TRUE},
then \code{"<link>(<name>)"} is returned instead of \code{"<name>"}.}
}
\value{
A character vector giving the subset of names relevant to \code{object},
or the complete set if \code{object = \link{NULL}}.
}
\description{
Retrieves the names used internally for top level nonlinear model
parameters.
}
