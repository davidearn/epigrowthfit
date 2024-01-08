\name{fixef.egf}
\alias{fixef.egf}
\alias{fixef}
\alias{fixef.egf_no_fit}
\title{Extract fixed effect coefficients}
\usage{
\method{fixef}{egf}(object, ...)

\method{fixef}{egf_no_fit}(object, ...)
}
\arguments{
\item{object}{An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.}

\item{...}{Unused optional arguments.}
}
\value{
A data frame with one row per coefficient and variables:
\item{bottom}{
  Name of a bottom level mixed effects model parameter,
  in this case a fixed effects coefficient;
  this is a string of the form \code{"beta[\%d]"}.
}
\item{top}{
  Name of the top level nonlinear model parameter whose fitted value
  is a function of \code{bottom},
  from \code{\link{egf_get_names_top}(object, link = TRUE)}.
}
\item{term}{
  Term from the fixed effects component of the mixed effects model formula
  for parameter \code{top}.
}
\item{colname}{
  Column name in the fixed effects design matrix
  \code{\link[=model.matrix.egf]{model.matrix}(object, "fixed")}.
}
\item{estimate}{
  Coefficient estimate, from segment \code{beta} of
  \code{\link[=coef.egf]{coef}(object, full = TRUE)}.
}
}
\description{
Retrieve the coefficients of the fixed effects component
of a mixed effects model.
}
\seealso{
Other extractors: 
\code{\link{coef.egf}()},
\code{\link{df.residual.egf}()},
\code{\link{fitted.egf}()},
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
