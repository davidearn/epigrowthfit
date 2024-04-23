\name{egf-class}
\alias{egf-class}
\title{Description of Objects of Class \code{egf}}
\description{
Class \code{egf} designates models estimated by function \code{\link{egf}}.
Objects of this class hold information about an estimated model.
Components can be accessed directly.  However, as the components
are subject to change without notice, portable code will rely on
exported methods for interrogation.
}
\details{
Currently, a legitimate \code{egf} object is a list with elements:

\describe{
\item{\code{model}}{
  a copy of the so-named argument of \code{\link{egf}}.}
\item{\code{frame}}{
  a list of the form \code{list(ts, windows, parameters, extra)}.
  \code{ts} and \code{windows} are data frames preserving time series
  and fitting window endpoints.
  \code{parameters} is a list of mixed effects model frames,
  with one element for each top level nonlinear model parameter.
  \code{extra} is a data frame preserving additional variables
  specified in \code{call[["select_windows"]]}.
  \code{windows}, the model frames listed in \code{parameters},
  and \code{extra} all correspond rowwise.}
\item{\code{priors}}{
  a list of the form \code{list(top, bottom = list(beta, theta, Sigma))},
  where \code{top}, \code{beta}, \code{theta}, and \code{Sigma} are all
  lists of \code{\link{egf_prior}} objects.}
\item{\code{control}}{
  a copy of the so-named argument of \code{\link{egf}}.}
\item{\code{tmb_out}}{
  the list output of \code{\link[TMB]{MakeADFun}}.}
\item{\code{optimizer_out}}{
  the list output of the optimizer specified by \code{control[["optimizer"]]}.}
\item{\code{init, best}}{
  numeric vectors giving the values of the condensed bottom level
  parameter vector used in the first and maximal likelihood evaluations.}
\item{\code{random}}{
  a logical vector indexing the elements of the condensed bottom level
  parameter vector that are \emph{not} arguments of the negative log
  marginal likelihood function.  It indexes
  all elements of segment \code{b} (random effect coefficients) and
  (but only if \code{control[["profile"]] = TRUE})
  all elements of segment \code{beta} (fixed effect coefficients).}
\item{\code{value, gradient}}{
  numeric vectors giving the value and gradient of the negative log
  marginal likelihood function at \code{best[!random]}.}
\item{\code{hessian}}{
  a logical flag indicating whether the Hessian matrix of the negative log
  marginal likelihood function is positive definite at \code{best[!random]}.
  \code{NA} means that the matrix has not been computed.}
\item{\code{coefficients}}{
  a list of the form \code{list(fixed, random)}, where \code{fixed}
  and \code{random} are data frames preserving interpretive information
  about fixed and random effect coefficients.}
\item{\code{contrasts}}{
  a list of the form \code{list(fixed, random)}, where \code{fixed}
  and \code{random} are lists preserving contrasts used to construct
  the fixed and random effects design matrices.}
\item{\code{call}}{
  the call to \code{\link{egf}}, enabling updates to the object by
  the default method of generic function \code{\link{update}}.}
}
}
\section{Bottom Level Parameter Vector}{
An estimated model is specified by a bottom level parameter vector
that is the concatenation of three segments:

\describe{
\item{\code{beta}}{
  the result of \code{unlist(lbeta)}, where \code{lbeta} is a list
  of numeric vectors of fixed effect coefficients,
  with one vector for each top level nonlinear model parameter.
  The order of top level parameters is specified by
  \code{\link{egf_top}(model)}.}
\item{\code{theta}}{
  the result of \code{unlist(ltheta)}, where \code{ltheta} is a list
  of numeric vectors of random effect covariance parameters,
  with one vector for each distinct random effect term in
  \code{formula_parameters}.
  Each vector parametrizes a random effect covariance matrix via
  \code{\link{theta2cov}} and its inverse \code{\link{cov2theta}}.\cr
  \cr
  The list \code{Sigma} mentioned in the description of \code{\link{egf}}
  argument \code{formula_priors} is precisely
  \code{\link{lapply}(ltheta, \link{theta2cov})}.}
\item{\code{b}}{
  the result of \code{unlist(lb)}, where \code{lb} is a list
  of numeric matrices of scaled random effect coefficients,
  corresponding elementwise to \code{ltheta}.
  The columns of \code{lb[[i]]} (one per level of the grouping variable)
  are interpreted as samples from a zero mean, unit variance multivariate
  normal distribution with covariance matrix
  \code{\link{cov2cor}(\link{theta2cov}(ltheta[[i]]))}.}
}

When elements of this vector are \dQuote{mapped} via \code{\link{egf}}
argument \code{map}, likelihood is defined as a function of the condensed
vector that excludes mapped elements.

Methods are defined for generic functions
\code{\link[=coef.egf]{coef}},
\code{\link[=fixef.egf]{fixef}}, and
\code{\link[=ranef.egf]{ranef}}
to allow users to interrogate the structure of the vector.
}
\examples{
methods(class = "egf")
help.search("\\\\.egf$", fields = "alias", package = "epigrowthfit")
## less verbosely: alias??`\\\\.egf$`
}
