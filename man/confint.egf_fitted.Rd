\name{confint.fitted.egf}
\alias{confint.fitted.egf}
\title{Confidence intervals on fitted values}
\usage{
\method{confint}{fitted.egf}(object, parm, level = 0.95, link = TRUE, ...)
}
\arguments{
\item{object}{An \code{"\link{fitted.egf}"} object.
Must supply link scale fitted values and corresponding standard errors.}

\item{parm}{Unused argument included for generic consistency.}

\item{level}{A number in the interval (0,1) indicating a confidence level.}

\item{link}{A \link{logical} flag. If \code{FALSE}, then confidence intervals
on inverse link-transformed fitted values are returned.}

\item{...}{Unused optional arguments.}
}
\value{
If \code{link = TRUE}, then \code{object} but with variable
\code{se} replaced with variables \code{lower} and \code{upper}
supplying confidence limits on fitted values (link scale).

Otherwise, the same result but with variables \code{estimate},
\code{lower}, and \code{upper} inverse link-transformed and
the \link{levels} of variable \code{top} modified accordingly.

\code{level} is retained as an \link[=attributes]{attribute}.
}
\description{
Computes confidence intervals on fitted values of top level nonlinear model
parameters.
}
\details{
Confidence limits on fitted values (link scale) are computed
as \code{estimate[i] + c(-1, 1) * sqrt(q) * se[i]},
with \code{estimate} and \code{se} as in \code{object} and
\code{q = \link{qchisq}(level, df = 1)}.
}
\examples{
object <- egf_cache("fitted-egf-1.rds")
zz <- confint(object)
str(zz)

}
