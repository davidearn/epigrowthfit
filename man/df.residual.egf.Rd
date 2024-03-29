\name{df.residual.egf}
\alias{df.residual.egf}
\alias{df.residual.egf_no_fit}
\title{Extract residual degrees of freedom}
\usage{
\method{df.residual}{egf}(object, ...)

\method{df.residual}{egf_no_fit}(object, ...)
}
\arguments{
\item{object}{An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.}

\item{...}{Unused optional arguments.}
}
\value{
An integer.
}
\description{
Returns the number of observations (see \code{\link[=nobs.egf]{nobs}})
minus the number of estimated parameters (fixed effects coefficients
and random effect covariance parameters).
}
\seealso{
Other extractors: 
\code{\link{coef.egf}()},
\code{\link{fitted.egf}()},
\code{\link{fixef.egf}()},
\code{\link{formula.egf}()},
\code{\link{getCall.egf}()},
\code{\link{logLik.egf}()},
\code{\link{model.frame.egf}()},
\code{\link{model.matrix.egf}()},
\code{\link{nobs.egf}()},
\code{\link{ranef.egf}()},
\code{\link{terms.egf}()},
\code{\link{vcov.egf}()}
}
\concept{extractors}
