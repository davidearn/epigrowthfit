\name{summary.egf}
\alias{summary.egf}
\title{Summarize a model object}
\usage{
\method{summary}{egf}(object, ...)
}
\arguments{
\item{object}{An \code{"\link{egf}"} object.}

\item{...}{Unused optional arguments.}
}
\value{
A list inheriting from class \code{"egf_summary"}, with elements:
\item{fitted}{
  A numeric matrix. Each column is the result of applying
  \code{\link{summary.default}} to a numeric vector listing
  the fitted values of a top level nonlinear model parameter.
  Fitted values are obtained via \code{\link[=fitted.egf]{fitted}}.
}
\item{convergence}{
  An integer code returned by the optimizer, with 0 indicating successful
  convergence within the specified absolute or relative tolerance.
}
\item{value, gradient}{
  Numeric vectors giving the value and gradient of the negative log
  likelihood function at the parameter vector returned by the optimizer.
}
\item{hessian}{
  A logical flag indicating whether the Hessian matrix of the negative log
  likelihood function at the parameter vector returned by the optimizer is
  positive definite.
  \code{NA} means that the matrix was not computed by \code{\link{egf}},
  either because \code{\link{egf}} was not called with \code{se = TRUE},
  or because an error was thrown during computation.
  In the first case, \code{object$sdreport} is \code{NULL}.
  In the second case, it is a \code{"\link[=try]{try-error}"} object
  preserving the error message.
}
}
\description{
Summarizes fitted values of top level nonlinear model parameters
and gathers diagnostic information that can be used to quickly
assess convergence of the optimizer.
}
\examples{
object <- egf_cache("egf-1.rds")
zz <- egf_cache("summary-egf-1.rds", summary(object))
str(zz)

}