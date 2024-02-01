\name{egf-class}
\alias{egf-class}
\title{Object of class "egf" or "egf_no_fit"}
\value{
A legitimate \code{"egf"} or \code{"egf_no_fit"} object is a list
with elements:
\item{model}{
  A copy of the so-named argument of \code{\link{egf}}.
}
\item{frame}{
  A list of the form \code{list(ts, windows, parameters, extra)}.
  \code{ts} and \code{windows} are data frames preserving time series
  and fitting window endpoints.
  \code{parameters} is a list of mixed effects model frames,
  with one element for each top level nonlinear model parameter.
  \code{extra} is a data frame preserving additional variables
  specified in \code{call$select_windows}.
  \code{windows}, the model frames listed in \code{parameters},
  and \code{extra} all correspond rowwise.
}
\item{priors}{
  A list of the form \code{list(top, bottom = list(beta, theta, Sigma))},
  where \code{top}, \code{beta}, \code{theta}, and \code{Sigma} are all
  lists of \code{"\link{egf_prior}"} objects.
}
\item{control}{
  A copy of the so-named argument of \code{\link{egf}}.
}
\item{tmb_out}{
  The list output of \code{\link[TMB]{MakeADFun}}.
  This contains an \link{environment} \code{env} whose objects
  are updated with each evaluation of the objective function.
}
\item{optimizer_out}{
  The list output of the optimizer specified by \code{control$optimizer}.
}
\item{init, best}{
  Numeric vectors giving the values of the condensed bottom level parameter
  vector used in the first and maximal likelihood evaluations.
  The \code{\link{names}} attribute of each vector groups the elements by
  segment; see Details.
}
\item{random}{
  A logical vector indexing the elements of the condensed bottom level
  parameter vector that are \emph{not} arguments of the negative log
  \emph{marginal} likelihood function.
  It indexes all elements of segment \code{b} (random effects) and,
  if \code{control$profile = TRUE}, all elements of segment \code{beta};
  see Details.
}
\item{value, gradient}{
  Numeric vectors giving the value and gradient of the negative log
  marginal likelihood function at \code{best[!random]}.
  \code{value} can be extracted using \code{\link[=logLik.egf]{logLik}}.
}
\item{hessian}{
  A logical flag indicating whether the Hessian matrix of the negative log
  (marginal) likelihood function is positive definite at \code{best[!random]}.
  \code{NA} means the matrix has not been computed
  (\code{sdreport} is not an \code{"sdreport"} object).
}
\item{sdreport}{
  If \code{call} contains \code{se = TRUE},
  then the result of \code{\link{try}(\link{sdreport}(tmb_out))}.
  Otherwise, \code{NULL}.
}
\item{effects}{
  A list of the form \code{list(beta, b)} with \code{beta} and \code{b}
  data frames preserving interpretive information about the so-named
  segments of the bottom level parameter vector.
}
\item{contrasts}{
  A list of the form \code{list(X, Z)}, with \code{X} and \code{Z}
  lists preserving the contrasts used to construct the fixed and
  random effects \link[=model.matrix]{design} matrices.
}
\item{call}{
  The \link{call} to \code{\link{egf}},
  enabling updates to the object via \code{\link{update}}.
}
\code{optimizer_out}, \code{best}, \code{value}, \code{gradient},
\code{hessian}, and \code{sdreport} have the value \code{NULL} in
objects of class \code{"egf_no_fit"}.
}
\description{
An object returned by \code{\link{egf}}, inheriting from class \code{"egf"}
or \code{"egf_no_fit"} depending on whether the nonlinear mixed effects model
specified in the call was actually fit.
}
\details{
Only developers should need to access the list directly,
e.g., using \code{$}.
Typical users can rely on methods to get information about the estimated
(or to-be-estimated) model.
See Examples for ways to query available methods and their documentation.

The estimated (or to-be-estimated) model is specified by a bottom level
parameter vector that is the concatenation of three segments:
\describe{
\item{beta}{
  The result of \code{unlist(lbeta)}, where \code{lbeta} is a list of
  numeric vectors of fixed effects coefficients, with one vector for
  each top level nonlinear model parameter.
  The order of top level parameters is given by
  \code{\link{egf_get_names_top}(model)}.
}
\item{theta}{
  The result of \code{unlist(ltheta)}, where \code{ltheta} is a list of
  numeric vectors of random effect covariance parameters, with one vector
  for each distinct random effects term in \code{formula_parameters}.
  Each vector parametrizes a random effect covariance matrix via
  \code{\link{theta2cov}} and its inverse \code{\link{cov2theta}}.
  (The list \code{Sigma} described under \code{\link{egf}} argument
  \code{formula_priors} is precisely
  \code{\link{lapply}(ltheta, \link{theta2cov})}.)
}
\item{b}{
  The result of \code{unlist(lb)}, where \code{lb} is a list of numeric
  matrices of scaled random effects coefficients, corresponding elementwise
  to \code{ltheta}.
  The columns of \code{lb[[i]]} (one per level of the grouping variable)
  are interpreted as samples from a zero mean, unit variance multivariate
  normal distribution with correlation matrix
  \code{\link{cov2cor}(\link{theta2cov}(ltheta[[i]]))}.
}
}
When elements of this vector are mapped via \code{\link{egf}} argument
\code{map}, likelihood is defined as a function of the condensed vector
that excludes mapped elements.

A number of methods are available to allow users to investigate the
structure of each of the three segments, with or without fitting a model;
see
\code{\link[=coef.egf]{coef}},
\code{\link[=fixef.egf]{fixef}}, and
\code{\link[=ranef.egf]{ranef}}.
}
\examples{
methods(class = "egf")
help.search("\\\\.egf$", fields = "alias", package = "epigrowthfit")
## less verbosely: alias??`\\\\.egf$`

methods(class = "egf_no_fit")
help.search("\\\\.egf_no_fit$", fields = "alias", package = "epigrowthfit")
## less verbosely: alias??`\\\\.egf_no_fit$`

}
