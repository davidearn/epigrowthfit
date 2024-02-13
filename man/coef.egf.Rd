\name{coef.egf}
\alias{coef.egf}
\alias{coef.egf_no_fit}
\title{Extract coefficients and covariance parameters}
\usage{
\method{coef}{egf}(object, full = FALSE, ...)

\method{coef}{egf_no_fit}(object, full = FALSE, ...)
}
\arguments{
\item{object}{An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.}

\item{full}{A logical flag. If \code{FALSE}, then the parameter vector is returned
in the condensed format used by \pkg{TMB}, which excludes mapped elements.
(See argument \code{map} of \code{\link{egf}}.)}

\item{...}{Unused optional arguments.}
}
\value{
A named numeric vector concatenating \code{beta}, \code{theta}, and \code{b},
without mapped elements if \code{full = FALSE}. The names group the elements
by segment.

Attribute \code{map} is a named list of index vectors \code{i} such that
that a full vector \code{y} and its condensed counterpart \code{x}
are related by \code{y = x[i]}, with the exception that \code{map[i]} is
\code{NA} if \code{y[i]} is mapped to an initial value
(hence \code{x[map[i]]} is equal to \code{NA}, not \code{y[i]}).
\code{NULL} is used in place of an index vector where \code{x} and \code{y}
are identical.

In both cases, the result inherits from class \code{"coef.egf"},
which has \code{\link{print}} and \code{\link{as.list}} methods.
}
\description{
Extracts the bottom level parameter vector \code{c(beta, theta, b)};
see \code{\link{egf}}.
}
\details{
If \code{object} inherits from class \code{"egf_no_fit"},
then the result contains initial parameter values to be used
in the first likelihood evaluation.
}
\seealso{
Other extractors: 
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
\code{\link{terms.egf}()},
\code{\link{vcov.egf}()}
}
\concept{extractors}
