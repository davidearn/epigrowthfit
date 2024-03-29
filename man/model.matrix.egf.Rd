\name{model.matrix.egf}
\alias{model.matrix.egf}
\alias{model.matrix.egf_no_fit}
\title{Extract design matrices}
\usage{
\method{model.matrix}{egf}(
  object,
  which = c("fixed", "random"),
  top = NULL,
  random = NULL,
  ...
)

\method{model.matrix}{egf_no_fit}(
  object,
  which = c("fixed", "random"),
  top = NULL,
  random = NULL,
  ...
)
}
\arguments{
\item{object}{An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.}

\item{which}{A character string controlling what is returned:
\code{"fixed"}, the fixed effects design matrix corresponding
to parameter \code{top}; or
\code{"random"}, the random effects design matrix corresponding
to parameter \code{top} and term \code{random}.}

\item{top}{A character string specifying a top level nonlinear model parameter.}

\item{random}{A random effects term, i.e., a call to the \code{|} operator.
Unused if \code{which = "fixed"}.}

\item{...}{Unused optional arguments.}
}
\value{
A \link[=matrix]{dense} or \link[Matrix:sparseMatrix]{sparse} matrix
with attributes \code{assign} and \code{contrasts}.
\code{assign} is absent in special cases; see Details.
}
\description{
Extract fixed and random effects \link[=model.matrix]{design} matrices
from a model object.
}
\details{
\code{model.matrix(which = "fixed", top = NULL)} returns the result of
combining (in the sense of \code{\link{cbind}}) all parameter-specific
fixed effects design matrices.

\code{model.matrix(which = "random", top = <string>, random = NULL)}
returns the result of combining (in the sense of \code{\link{cbind}})
all term-specific random effects design matrices associated with
parameter \code{top}.

\code{model.matrix(which = "random", top = NULL)} returns the result
of combining the parameter-specific all-terms matrices \emph{and}
permuting the columns to obtain a convenient ordering of random effects
coefficients.
(Coefficients are sorted by relation to a common random vector;
random vectors are sorted by relation to a common covariance matrix.)

None of these combined matrices are "proper" design matrices,
as none possess an \code{assign} attribute.
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
\code{\link{nobs.egf}()},
\code{\link{ranef.egf}()},
\code{\link{terms.egf}()},
\code{\link{vcov.egf}()}
}
\concept{extractors}
