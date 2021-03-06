% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitted.R
\name{fitted.egf}
\alias{fitted.egf}
\alias{coef.egf}
\title{Extract fitted values}
\usage{
\method{fitted}{egf}(
  object,
  par = get_par_names(object, link = TRUE),
  subset = NULL,
  append = NULL,
  link = TRUE,
  se = FALSE,
  .subset = NULL,
  .append = NULL,
  ...
)

\method{coef}{egf}(
  object,
  par = get_par_names(object, link = TRUE),
  subset = NULL,
  append = NULL,
  link = TRUE,
  se = FALSE,
  .subset = NULL,
  .append = NULL,
  ...
)
}
\arguments{
\item{object}{An \code{"\link{egf}"} object.}

\item{par}{A subset of \code{\link{get_par_names}(object, link = TRUE)}
naming nonlinear and dispersion model parameters for which
fitted values should be retrieved.}

\item{subset}{An expression to be evaluated in the combined model frame
(see \code{\link{make_combined}}). Must evaluate to
a \link{logical} vector indexing rows of the data frame,
and thus fitting windows. Fitted values are retrieved
only for indexed windows. The default (\code{\link{NULL}})
is to consider all windows.}

\item{append}{An expression indicating variables in the combined model frame
(see \code{\link{make_combined}}) to be included with the result.
The default (\code{\link{NULL}}) is to append nothing.}

\item{link}{A \link{logical} flag. If \code{FALSE}, then fitted values
are inverse link-transformed.}

\item{se}{A \link{logical} flag. If \code{link = TRUE} and \code{se = TRUE},
then approximate (delta method) standard errors on fitted values
are reported. Standard errors are required for subsequent use of
\code{\link{confint.egf_fitted}}.}

\item{.subset}{A \link{logical} vector to be used (if non-\code{\link{NULL}})
in place of the result of evaluating \code{subset}.}

\item{.append}{A \link{character} vector listing variable names to be used
(if non-\code{\link{NULL}}) in place of the result of evaluating
\code{append}.}

\item{...}{Unused optional arguments.}
}
\value{
A \link[=data.frame]{data frame} inheriting from \link{class}
\code{"egf_fitted"}, with variables:
\item{par}{
  Nonlinear or dispersion model parameter,
  from \code{\link{get_par_names}(object, link = link)}.
}
\item{ts}{
  Time series, from \code{\link{levels}(object$endpoints$ts)}.
}
\item{window}{
  Fitting window, from \code{\link{levels}(object$endpoints$window)}.
}
\item{estimate}{
  Fitted value of parameter \code{par} in fitting window \code{window}.
}
\item{se}{
  (If \code{link = TRUE} and \code{se = TRUE}.)
  Approximate (delta method) standard error on \code{estimate}.
}
}
\description{
Computes fitted values of nonlinear and dispersion model parameters.
The fitted value for a given fitting window is obtained by adding
(i) the population fitted value computed as a linear combination
of fixed effects coefficients and
(ii) all applicable random effects, with random effects coefficients
set equal to their conditional modes.
}
\details{
\code{coef.egf} is currently an alias for \code{fitted.egf}.

See topic \code{\link{nse}} for details on nonstandard evaluation
of \code{subset} and \code{append}.
}
\seealso{
\code{\link{confint.egf_fitted}}
}
