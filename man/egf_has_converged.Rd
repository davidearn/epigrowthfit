\name{egf_has_converged}
\alias{egf_has_converged}
\title{Check for convergence}
\usage{
egf_has_converged(object, tol = 1)
}
\arguments{
\item{object}{An \code{"\link{egf}"} object.}

\item{tol}{A positive number. Convergence requires all gradient elements
to be less than or equal to \code{tol} in absolute value.}
}
\value{
\code{TRUE} if all tests pass. \code{FALSE} if any test fails.
\code{NA} if no test fails, but the test for a positive definite
Hessian matrix is indeterminate because the matrix has not been
computed.
}
\description{
Performs simple diagnostic checks to assess whether the optimizer that
produced an estimated model actually converged to a local minimum point
of the negative log likelihood function.
}
