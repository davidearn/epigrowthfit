\name{logLik.egf}
\alias{logLik.egf}
\title{Extract log marginal likelihood}
\usage{
\method{logLik}{egf}(object, ...)
}
\arguments{
\item{object}{An \code{"\link{egf}"} object.}

\item{...}{Unused optional arguments.}
}
\value{
A number inheriting from class \code{"logLik"}.
Attribute \code{df} is the number of estimated parameters
(fixed effects coefficients and random effect covariance parameters).
Attribute \code{nobs} is the number of observations of disease incidence
used in estimation.
}
\description{
Retrieves the value of log marginal likelihood from a fitted model object.
}
\seealso{
Other extractors: 
\code{\link{coef.egf}()},
\code{\link{df.residual.egf}()},
\code{\link{fitted.egf}()},
\code{\link{fixef.egf}()},
\code{\link{formula.egf}()},
\code{\link{getCall.egf}()},
\code{\link{model.frame.egf}()},
\code{\link{model.matrix.egf}()},
\code{\link{nobs.egf}()},
\code{\link{ranef.egf}()},
\code{\link{terms.egf}()},
\code{\link{vcov.egf}()}
}
\concept{extractors}
