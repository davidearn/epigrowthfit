\name{terms.egf}
\alias{terms.egf}
\alias{terms.egf_no_fit}
\title{Extract model terms}
\usage{
\method{terms}{egf}(x, top = egf_get_names_top(x, link = TRUE), ...)

\method{terms}{egf_no_fit}(x, top = egf_get_names_top(x, link = TRUE), ...)
}
\arguments{
\item{x}{An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.}

\item{top}{A character string specifying a top level nonlinear model parameter.}

\item{...}{Unused optional arguments.}
}
\value{
A \code{"\link[=terms.object]{terms}"} object.
}
\description{
Extract \code{"\link[=terms.object]{terms}"} objects
corresponding to top level nonlinear model parameters.
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
\code{\link{vcov.egf}()}
}
\concept{extractors}
