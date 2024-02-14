\name{summary.egf}
\alias{summary.egf}
\title{Summarize a Model Object}
\description{
Summarizes fitted values of top level nonlinear model parameters
and gathers diagnostic information that can be used to quickly assess
convergence of the optimizer.
}
\usage{
\method{summary}{egf}(object, \dots)
}
\arguments{
\item{object}{an \code{\link{egf}} object.}
\item{\dots}{unused optional arguments.}
}
\value{
A list inheriting from class \code{summary.egf}, with elements:

\item{fitted}{a numeric matrix.  Each column is the result of applying
  \code{\link{summary.default}} to a numeric vector listing the fitted
  values of a top level nonlinear model parameters.  Fitted values are
  retrieved by \code{\link{fitted.egf}}.}
\item{convergence}{an integer code returned by the optimizer, with 0
  indicating successful convergence within the specified absolute or
  relative tolerance.}
\item{value, gradient}{numeric vectors giving the value and gradient
  of the negative log likelihood function at the parameter vector
  returned by the optimizer.}
\item{hessian}{a logical flag indicating whether the Hessian matrix
  of the negative log likelihood function is positive definite at the
  parameter vector returned by the optimizer.  \code{NA} means that
  the matrix was not computed by \code{\link{egf}}, either because
  \code{se = TRUE} was not passed in the function call or because an
  error was thrown during computation.}
}
\examples{
object <- egf_cache("egf-1.rds")
zz <- egf_cache("summary-egf-1.rds", summary(object))
str(zz)
}
