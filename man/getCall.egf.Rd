\name{getCall.egf}
\alias{getCall.egf}
\alias{getCall.egf_no_fit}
\title{Extract model calls}
\usage{
\method{getCall}{egf}(x, ...)

\method{getCall}{egf_no_fit}(x, ...)
}
\arguments{
\item{x}{An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.}

\item{...}{Unused optional arguments.}
}
\value{
A \link{call} to \code{\link{egf}}.
}
\description{
Extract the \link{call} to \code{\link{egf}} that produced the given
model object. These methods exist mainly to enable compatibility with
the default method of \code{\link{update}}.
}
\seealso{
Other extractors: 
\code{\link{coef.egf}()},
\code{\link{df.residual.egf}()},
\code{\link{fitted.egf}()},
\code{\link{fixef.egf}()},
\code{\link{formula.egf}()},
\code{\link{logLik.egf}()},
\code{\link{model.frame.egf}()},
\code{\link{model.matrix.egf}()},
\code{\link{nobs.egf}()},
\code{\link{ranef.egf}()},
\code{\link{terms.egf}()},
\code{\link{vcov.egf}()}
}
\concept{extractors}
