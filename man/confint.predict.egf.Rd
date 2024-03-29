\name{confint.predict.egf}
\alias{confint.predict.egf}
\title{Confidence intervals on predicted values}
\usage{
\method{confint}{predict.egf}(object, parm, level = 0.95, log = TRUE, ...)
}
\arguments{
\item{object}{An \code{"\link{predict.egf}"} object.
Must supply log scale predicted values and corresponding standard errors.}

\item{parm}{Unused argument included for generic consistency.}

\item{level}{A number in the interval (0,1) indicating a confidence level.}

\item{log}{A logical flag. If \code{FALSE}, then confidence intervals
on inverse log-transformed predicted values are returned.}

\item{...}{Unused optional arguments.}
}
\value{
If \code{log = TRUE}, then \code{object} but with variable
\code{se} replaced with variables \code{lower} and \code{upper}
supplying confidence limits on log predicted values.

Otherwise, the same result but with variables \code{estimate},
\code{lower}, and \code{upper} inverse log-transformed and
the \link{levels} of variable \code{var} modified accordingly.

\code{level} is retained as an attribute.
}
\description{
Computes confidence intervals on predicted values of interval incidence,
cumulative incidence, and the per capita growth rate.
}
\details{
Confidence limits on predicted values (log scale) are computed
as \code{estimate + c(-1, 1) * sqrt(q) * se},
with \code{estimate} and \code{se} as in \code{object} and
\code{q = \link{qchisq}(level, df = 1)}.
}
\examples{
object <- egf_cache("predict-egf-1.rds")
confint(object, log = TRUE)
confint(object, log = FALSE)

}
