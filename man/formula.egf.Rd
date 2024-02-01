\name{formula.egf}
\alias{formula.egf}
\alias{formula.egf_no_fit}
\title{Extract model formulae}
\usage{
\method{formula}{egf}(x, top = egf_top(x, link = TRUE), split = FALSE, ...)

\method{formula}{egf_no_fit}(x, top = egf_top(x, link = TRUE), split = FALSE, ...)
}
\arguments{
\item{x}{An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.}

\item{top}{A character string specifying a top level nonlinear model parameter.}

\item{split}{A logical flag. If \code{TRUE}, then fixed and random effects terms
are returned separately.}

\item{...}{Unused optional arguments.}
}
\value{
If \code{split = FALSE}, then the mixed effects model formula
corresponding to parameter \code{top}.
If \code{split = TRUE}, then the same formula but with random effects
terms deleted from the right hand side.
The deleted terms are preserved in a list and retained as an attribute
of the result, namely \code{random}.
}
\description{
Extract mixed effects model formulae corresponding to top level
nonlinear model parameters.
}
\seealso{
Other extractors: 
\code{\link{coef.egf}()},
\code{\link{df.residual.egf}()},
\code{\link{fitted.egf}()},
\code{\link{fixef.egf}()},
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
