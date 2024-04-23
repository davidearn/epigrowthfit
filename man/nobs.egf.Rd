\name{nobs.egf}
\alias{nobs.egf}
\alias{nobs.egf_no_fit}
\title{Extract the Number of Observations}
\description{
Returns the number of observations of disease incidence that were
used in estimation of a model.  This number excludes missing values
and observations not belonging to a fitting window, which,
despite being preserved in model objects, do not affect estimation.
}
\usage{
\method{nobs}{egf}(object, \dots)
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
The generic function \code{\link{nobs}}.
}
