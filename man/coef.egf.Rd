\name{coef.egf}
\alias{coef.egf}
\alias{coef.egf_no_fit}
\alias{print.coef.egf}
\alias{as.list.coef.egf}
\title{Extract Coefficients and Random Effect Covariance Parameters}
\description{
Extracts the bottom level parameter vector \code{c(beta, theta, b)}
or a subset.  Segments \code{beta}, \code{theta}, and \code{b} contain
(respectively) fixed effect coefficients, random effect covariance
parameters, and random effect coefficients.
}
\usage{
\method{coef}{egf}(object, random = FALSE, full = FALSE, \dots)
\method{coef}{egf_no_fit}(object, \dots)
}
\arguments{
\item{object}{an \code{\link{egf}} or \code{\link[=egf]{egf_no_fit}}
  object.}
\item{random}{a logical.  If \code{FALSE}, then segment \code{b}
  is excluded.}
\item{full}{a logical.  If \code{FALSE}, then mapped elements are
  excluded, and the result is called \dQuote{condensed}.}
\item{\dots}{unused optional arguments.}
}
\value{
A named numeric vector concatenating \code{beta}, \code{theta}, and
\code{b}, without \code{b} if \code{random = FALSE} and without
mapped elements if \code{full = FALSE}.  The names group the elements
by segment.

Attribute \code{map} is a named list of integer vectors \code{i}
such that that a full segment \code{y} and its condensed counterpart
\code{x} are related by \code{y = x[i]}, with the exception that
\code{i[j]} is \code{NA} if \code{y[j]} is mapped to an initial value.
\code{NULL} is used in place of an integer vector where \code{x}
and \code{y} are identical.

The result inherits from class \code{coef.egf}, which has methods
for \code{\link{print}} and \code{\link{as.list}}.
}
\note{
The method for class \code{egf_no_fit} returns the initial value
of the parameter vector, not an estimated value.
}
\examples{
object <- egf_cache("egf-2.rds")
for (random in c(FALSE, TRUE)) {
    for (full in c(FALSE, TRUE)) {
        cat(sprintf("random = \%s, full = \%s :\n\n", random, full))
        str(coef(object, random = random, full = full))
        cat("\n")
    }
}
}
