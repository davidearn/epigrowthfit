% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update.R
\name{getCall.egf}
\alias{getCall.egf}
\alias{getCall.egf_no_fit}
\alias{getCall.egf_model_simulate}
\title{Update an object returned by an S3 method}
\usage{
\method{getCall}{egf}(x, ...)

\method{getCall}{egf_no_fit}(x, ...)

\method{getCall}{egf_model_simulate}(x, ...)
}
\arguments{
\item{x}{A \link{list} with an element \code{call}, which should be a \link{call}.}

\item{...}{Unused optional arguments.}
}
\value{
\code{x$call}, modified so that \code{x$call[[1L]]} is a generic function,
not an S3 method.
}
\description{
When \code{\link{match.call}} is called from an S3 method for a generic
function, a \link{call} to the method is returned, rather than a call
to the generic function. Typically, the method is not exported, so the
default method for \code{\link{update}} is unable to evaluate the call.
These methods for \code{\link{getCall}} circumvent this issue for objects
created by various S3 methods implemented in \pkg{epigrowthfit}.
}