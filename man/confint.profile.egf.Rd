\name{confint.profile.egf}
\alias{confint.profile.egf}
\title{Confidence intervals from likelihood profiles}
\usage{
\method{confint}{profile.egf}(object, parm, level = attr(object, "level"), link = TRUE, ...)
}
\arguments{
\item{object}{An \code{"\link{profile.egf}"} object.}

\item{parm}{Unused argument included for generic consistency.}

\item{level}{A number in the interval (0,1) indicating a confidence level.}

\item{link}{A logical flag. If \code{FALSE} and \code{object} supplies
likelihood profiles of population fitted values of top level
nonlinear model parameters, then confidence intervals
on inverse link-transformed fitted values are returned.}

\item{...}{Unused optional arguments.}
}
\value{
A data frame with one row per supplied profile, and variables:
\item{linear_combination}{
  Row index of linear combination that was profiled,
  from \code{seq_len(nrow(attr(object, "A")))}.
}
\item{estimate, lower, upper}{
  Estimate of linear combination and approximate lower
  and upper confidence limits, inverse-link transformed
  if \code{link = FALSE}.
}
\code{level} is retained as an attribute of the result.
So are attributes \code{A} and \code{x} of \code{object}.
}
\description{
Computes confidence intervals on fixed effects coefficients,
random effects covariance parameters, and linear combinations
thereof from their univariate likelihood profiles.
}
\details{
Each supplied likelihood profile
(level of factor \code{object$linear_combination}),
is linearly interpolated to approximate the two solutions
of \code{deviance(value) = \link{qchisq}(level, df = 1)}.
These provide the lower and upper confidence limits of interest
(see \href{https://en.wikipedia.org/wiki/Wilks'_theorem}{Wilks' theorem}).
}
\examples{
object <- egf_cache("profile-egf-1.rds")
zz <- confint(object)
str(zz)

}
