\name{vcov.egf}
\alias{vcov.egf}
\title{Model Covariance Matrices}
\description{
Extracts (or, if necessary, computes) the covariance matrix of
bottom level parameters \code{beta} and \code{theta},
corresponding to the output of \code{\link[=coef.egf]{coef}(object)}.
}
\usage{
\method{vcov}{egf}(object, \dots)
}
\arguments{
\item{object}{
  an \code{\link{egf}} object.}
\item{\dots}{
  unused optional arguments.}
}
\value{
A real, symmetric matrix.
}
\details{
If the resulting matrix is not finite or not positive definite,
then the fit specified by \code{object} should be investigated,
as the optimizer that produced the fit may have failed to converge.
}
\seealso{
The generic function \code{\link{vcov}}.
}
