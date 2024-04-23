\name{logLik.egf}
\alias{logLik.egf}
\title{Extract the Log (Marginal) Likelihood}
\description{
Extracts from a model object the value of the log marginal likelihood.
Whether the result represents a local maximum depends on, among other
things, convergence of the optimizer.
}
\usage{
\method{logLik}{egf}(object, \dots)
}
\arguments{
\item{object}{
  an \code{\link{egf}} object.}
\item{\dots}{
  unused optional arguments.}
}
\value{
A numeric vector of length 1 inheriting from class \code{logLik}.
Attribute \code{df} is the number of estimated parameters
(fixed effect coefficients and random effect covariance parameters).
Attribute \code{nobs} is the number of observations of disease
incidence used in estimation.
}
\seealso{
The generic function \code{\link{logLik}}.
}
