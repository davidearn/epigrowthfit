\name{vcov.egf}
\alias{vcov.egf}
\title{Extract model covariance matrix}
\usage{
\method{vcov}{egf}(object, ...)
}
\arguments{
\item{object}{An \code{"\link{egf}"} object.}

\item{...}{Unused optional arguments.}
}
\value{
A symmetric matrix.
}
\description{
Extracts (or, if necessary, computes) the covariance matrix
of bottom level parameters \code{beta} and \code{theta},
corresponding to the output of
\code{\link[=coef.egf]{coef}(object, full = FALSE)}.
}
\details{
If \code{object} was constructed by a call \code{\link{egf}(se = TRUE)},
then the full covariance matrix has already been computed and is preserved
in \code{object}. \code{vcov} reuses this matrix to avoid recomputation.

If the returned matrix is not finite and positive definite,
then the fit specified by \code{object} should be investigated,
as the optimizer that produced the fit may have failed to converge.
See also \code{\link{egf_has_converged}}.
}
\seealso{
Other extractors: 
\code{\link{coef.egf}()},
\code{\link{df.residual.egf}()},
\code{\link{fitted.egf}()},
\code{\link{fixef.egf}()},
\code{\link{formula.egf}()},
\code{\link{getCall.egf}()},
\code{\link{logLik.egf}()},
\code{\link{model.frame.egf}()},
\code{\link{model.matrix.egf}()},
\code{\link{nobs.egf}()},
\code{\link{ranef.egf}()},
\code{\link{terms.egf}()}
}
\concept{extractors}
