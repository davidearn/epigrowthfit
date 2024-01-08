\name{profile.egf}
\alias{profile.egf}
\title{Compute likelihood profiles}
\usage{
\method{profile}{egf}(
  fitted,
  level = 0.95,
  which = NULL,
  A = NULL,
  top = egf_get_names_top(fitted, link = TRUE),
  parallel = egf_parallel(),
  trace = FALSE,
  grid_len = 12,
  subset = NULL,
  append = NULL,
  .subset = NULL,
  .append = NULL,
  ...
)
}
\arguments{
\item{fitted}{An \code{"\link{egf}"} object.}

\item{level}{A number in the interval (0,1) indicating a confidence level.
Profiles will be computed up to a deviance
of \code{\link{qchisq}(level, df = 1)}.}

\item{which}{An integer vector indexing \code{fitted$best[!fitted$random]}.
Only indexed parameters are profiled.}

\item{A}{A numeric matrix with \code{sum(!fitted$random)} columns,
Each row specifies a linear combination of the elements
of \code{fitted$best[!fitted$random]} to be profiled.
Ignored if \code{which} is non-\code{NULL}.}

\item{top}{A subset of \code{\link{egf_get_names_top}(fitted, link = TRUE)}
naming top level nonlinear model parameters whose population
fitted values (see \code{\link[=fitted.egf]{fitted}}) should be
profiled.
Ignored if \code{which} or \code{A} is non-\code{NULL}.}

\item{parallel}{An \code{"\link{egf_parallel}"} object defining options for \R level
parallelization.}

\item{trace}{A logical flag.
If \code{TRUE}, then basic tracing messages indicating progress
are printed.
Depending on \code{fitted$control$trace}, these may be mixed with
optimizer output.}

\item{grid_len}{A positive integer. Step sizes chosen adaptively by
\code{\link[TMB]{tmbprofile}} will generate approximately
this many points on each side of a profile's minimum point.}

\item{subset}{An expression to be evaluated in
\code{\link[=model.frame.egf]{model.frame}(fitted, "combined")}
It must evaluate to to a valid index vector for the rows of the
data frame and, in turn, fitting windows.
Only population fitted values for indexed windows are profiled.
The default (\code{NULL}) is to consider all windows.
Ignored if \code{which} or \code{A} is non-\code{NULL}.}

\item{append}{An expression indicating variables in
\code{\link[=model.frame.egf]{model.frame}(object, "combined")}
to be included with the result.
The default (\code{NULL}) is to append nothing.}

\item{.subset, .append}{Index vectors to be used (if non-\code{\link{NULL}}) in place of
the result of evaluating \code{subset} and \code{append}.}

\item{...}{Unused optional arguments.}
}
\value{
A data frame inheriting from class \code{"egf_profile"}, with variables:
\item{top}{
  (\code{top}-based calls only.)
  Top level nonlinear model parameter,
  from \code{\link{egf_get_names_top}(fitted, link = link)}.
}
\item{ts}{
  (\code{top}-based calls only.)
  Time series, from
  \code{levels(\link[=model.frame.egf]{model.frame}(fitted)$ts)}.
}
\item{window}{
  (\code{top}-based calls only.)
  Fitting window, from
  \code{levels(\link[=model.frame.egf]{model.frame}(fitted)$window)}.
}
\item{linear_combination}{
  Row index of linear combination, from \code{seq_len(nrow(A))}.
}
\item{value}{
  Value of linear combination being profiled.
}
\item{deviance}{
  Deviance of the restricted model that assumes \code{value}
  for the linear combination being profiled.
}
\code{A}, \code{x = fitted$best[!fitted$random]}, and \code{level}
are retained as \link{attributes}.
}
\description{
Computes univariate likelihood profiles of fixed effects
coefficients, random effect covariance parameters, and
linear combinations thereof.
}
\details{
Computation of likelihood profiles is expensive as it requires estimation
of many restricted models. It is parallelized at the C++ level when there
is OpenMP support and \code{fitted$control$omp_num_threads} is set to
an integer greater than 1. If there is no OpenMP support, then computation
can still be parallelized at the \R level with appropriate setting of
\code{parallel}.

See topic \code{\link{egf_eval}} for details on nonstandard evaluation
of \code{subset} and \code{append}.
}
\examples{
fitted <- egf_cache("egf-1.rds")
zz <- egf_cache("profile-egf-1.rds", {
  profile(fitted, subset = (country == "A" & wave == 1))
})
str(zz)

}
\seealso{
\code{\link{confint.egf_profile}}, \code{\link{plot.egf_profile}}
}
