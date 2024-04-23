\name{egf_has_converged}
\alias{egf_has_converged}
\title{Test for Convergence}
\description{
Performs simple diagnostic tests to assess whether the optimizer that
produced an estimated model actually converged to a local minimum point
of the negative log marginal likelihood function.
}
\usage{
egf_has_converged(object, check = TRUE, tol = 1)
}
\arguments{
\item{object}{
  an \code{\link{egf}} object.}
\item{check}{
  a logical.  If \code{TRUE}, then an error is thrown if \code{object}
  does not actually inherit from class \code{\link{egf}}.}
\item{tol}{
  a positive number.  Convergence requires all gradient elements
  to be less than or equal to \code{tol} in absolute value.}
}
\value{
\code{TRUE} if all tests pass.  \code{FALSE} if any test fails.
\code{NA} if no test fails, but the test for a positive definite
Hessian matrix is indeterminate because the matrix has not been
computed.
}
