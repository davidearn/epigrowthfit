\name{df.residual.egf}
\alias{df.residual.egf}
\alias{df.residual.egf_no_fit}
\title{Extract the Residual Degrees of Freedom}
\description{
Extracts from a model object the number of observations
(see \code{\link[=nobs.egf]{nobs}})
minus the number of estimated parameters
(fixed effect coefficients and random effect covariance parameters).
}
\usage{
\method{df.residual}{egf}(object, \dots)
}
\arguments{
\item{object}{
  an \code{\link{egf}} object.}
\item{\dots}{
  unused optional arguments.}
}
\value{
An integer.
}
\seealso{
The generic function \code{\link{df.residual}}.
}
