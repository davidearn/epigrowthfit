\name{nobs.egf}
\alias{nobs.egf}
\alias{nobs.egf_no_fit}
\title{Extract number of observations}
\usage{
\method{nobs}{egf}(object, ...)

\method{nobs}{egf_no_fit}(object, ...)
}
\arguments{
\item{object}{An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.}

\item{...}{Unused optional arguments.}
}
\value{
An integer.
}
\description{
Returns the number of observations of disease incidence that were used
or would be used in estimation of a model. This number excludes missing
values and observations not belonging to a fitting window, which,
despite being preserved in model objects, do not affect estimation.
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
\code{\link{ranef.egf}()},
\code{\link{terms.egf}()},
\code{\link{vcov.egf}()}
}
\concept{extractors}
