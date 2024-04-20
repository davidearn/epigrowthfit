\name{ranef.egf}
\alias{ranef.egf}
\alias{ranef}
\alias{ranef.egf_no_fit}
\title{Extract random effect conditional modes}
\usage{
\method{ranef}{egf}(object, makeSigma = FALSE, ...)

\method{ranef}{egf_no_fit}(object, makeSigma = FALSE, ...)
}
\arguments{
\item{object}{An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.}

\item{makeSigma}{A logical flag.
If \code{TRUE}, then random effect covariance matrices are constructed
from segment \code{theta} of
\code{\link[=coef.egf]{coef}(object, full = TRUE)}
and preserved as an attribute of the result.}

\item{...}{Unused optional arguments.}
}
\value{
A data frame with one row per coefficient and variables:
\item{cov}{
  Name of a covariance matrix.
  This is the interaction of \code{term} and \code{group},
  but with levels named \code{"Sigma[\%d]"}.
}
\item{vec}{
  Name of a random vector.
  This is the interaction of \code{term}, \code{group}, and \code{level},
  but with levels named \code{"u[\%d]"}.
}
\item{bottom}{
  Name of a bottom level mixed effects model parameter,
  in this case a random effects coefficient;
  this is a string of the form \code{"b[\%d]"}.
}
\item{top}{
  Name of the top level nonlinear model parameter whose fitted value
  is a function of \code{bottom}.
}
\item{term, group}{
  Random effects term from mixed effects model formula
  for parameter \code{top}. \code{term} and \code{group} give
  the left and right hand sides of the \code{`|`} operator.
}
\item{level}{
  Level of factor or interaction indicated by \code{group}.
}
\item{colname}{
  Column name in the random effects design matrix
  \code{\link[=model.matrix.egf]{model.matrix}(object, "random")}.
}
\item{mode}{
  Random effect conditional mode (unit variance scale),
  from segment \code{b} of
  \code{\link[=coef.egf]{coef}(object, full = TRUE)}.
}
If \code{makeSigma = TRUE}, then the result has attribute \code{Sigma},
a list of covariance matrices corresponding to the levels of variable
\code{cov}.
}
\description{
Retrieve the coefficients of the random effects component
of a mixed effects model
(specifically, their modes conditional on data and parameter estimates).
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
\code{\link{terms.egf}()},
\code{\link{vcov.egf}()}
}
\concept{extractors}
