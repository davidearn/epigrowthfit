\name{confint.egf}
\alias{confint.egf}
\title{Confidence intervals on fitted values}
\usage{
\method{confint}{egf}(
  object,
  parm,
  level = 0.95,
  top = egf_top(object),
  link = TRUE,
  method = c("wald", "profile", "uniroot"),
  parallel = egf_parallel(),
  trace = FALSE,
  grid_len = 12,
  interval_scale = 7,
  subset = NULL,
  select = NULL,
  ...
)
}
\arguments{
\item{object}{An \code{"\link{egf}"} object.}

\item{parm}{Unused argument included for generic consistency.}

\item{level}{A number in the interval (0,1) indicating a confidence level.}

\item{top}{A subset of \code{\link{egf_top}(object, link = TRUE)}
naming top level nonlinear model parameters for which confidence
intervals should be computed.}

\item{link}{A logical flag. If \code{FALSE}, then confidence intervals
on inverse link-transformed fitted values are computed.}

\item{method}{A \link{character} string indicating how confidence intervals
should be computed.}

\item{parallel}{(For \code{method = "profile"} or \code{"uniroot"}.)
An \code{"\link{egf_parallel}"} object defining options for \R level
parallelization.}

\item{trace}{(For \code{method = "profile"} or \code{"uniroot"}.)
A logical flag.
If \code{TRUE}, then basic tracing messages indicating progress
are printed.
Depending on \code{object$control$trace}, these may be mixed with
optimizer output.}

\item{grid_len}{(For \code{method = "profile"}.)
A positive integer. Step sizes chosen adaptively by
\code{\link[TMB]{tmbprofile}} will generate approximately
this many points on each side of a profile's minimum point.}

\item{interval_scale}{(For \code{method = "uniroot"}.)
A positive number.
\code{\link[TMB]{tmbroot}} will search for roots in the interval
of width \code{interval_scale * se} centered at \code{estimate},
where \code{estimate} is the fitted value (link scale)
and \code{se} is the corresponding standard error.}

\item{subset, select}{Index vectors for the rows and columns of
\code{mf = \link[=model.frame.egf]{model.frame}(object, "combined")}
or language objects evaluating to such vectors.
\code{subset} indicates fitting windows for which fitted values
should be retrieved; the default (\code{NULL}) indicates all of them.
\code{select} indicates variables that should be appended to the result;
the default (\code{NULL}) is to append nothing.
Evaluation of language objects follows the implementation of
\code{\link{subset}(mf, subset, select)}.}

\item{...}{Unused optional arguments.}
}
\value{
A data frame inheriting from class \code{"confint.egf"}, with variables:
\item{top}{
  Top level nonlinear model parameter,
  from \code{\link{egf_top}(object, link = TRUE)}.
}
\item{ts}{
  Time series, from
  \code{\link{levels}(\link[=model.frame.egf]{model.frame}(object)$ts)}.
}
\item{window}{
  Fitting window, from
  \code{\link{levels}(\link[=model.frame.egf]{model.frame}(object)$window)}.
}
\item{estimate, lower, upper}{
  Fitted value and approximate lower and upper confidence limits.
}
\code{level} and
\code{frame_windows = \link[=model.frame.egf]{model.frame}(object, "windows")}
are retained as attributes.
}
\description{
Computes confidence intervals on fitted values of top level nonlinear
model parameters.
}
\details{
Three methods are provided for calculating confidence intervals:
\describe{
\item{\code{wald}}{
  See \code{\link{confint.fitted.egf}}.
}
\item{\code{profile}}{
  See \code{\link{confint.profile.egf}}.
}
\item{\code{uniroot}}{
  Similar to \code{"profile"}, except that the two solutions
  of \code{deviance(value) = \link{qchisq}(level, df = 1)} are
  approximated by root-finding using \code{\link[TMB]{tmbroot}}
  (\code{\link{uniroot}} internally).
}
}
For top level parameters following random effects models,
\code{"wald"} returns confidence intervals on individual fitted values,
while \code{"profile"} and \code{"uniroot"} return confidence intervals
on population fitted values, which are the fixed effects components of
individual fitted values.

\code{"wald"} assumes, e.g., asymptotic normality of the maximum likelihood
estimator. \code{"profile"} and \code{"uniroot"} avoid these issues but are
typically more expensive, requiring estimation of many restricted models.
They are parallelized at the C++ level when there is OpenMP support and
\code{object$control$omp_num_threads} is set to an integer greater than 1.
If there is no OpenMP support, then computation can still be parallelized
at the \R level with appropriate setting of \code{parallel}.
}
\examples{
object <- egf_cache("egf-1.rds")
zz1 <- egf_cache("confint-egf-1.rds", {
    confint(object, method = "wald")
})
zz2 <- egf_cache("confint-egf-2.rds", {
    confint(object, method = "profile", top = "log(r)",
            subset = quote(country == "A" & wave == 1))
})
zz3 <- egf_cache("confint-egf-3.rds", {
    confint(object, method = "uniroot", top = "log(r)",
            subset = quote(country == "A" & wave == 1))
})
str(zz1)

}
\seealso{
\code{\link{plot.confint.egf}}
}
