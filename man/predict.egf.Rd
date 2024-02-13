\name{predict.egf}
\alias{predict.egf}
\title{Compute predicted values}
\usage{
\method{predict}{egf}(
  object,
  what = c("interval", "cumulative", "rt"),
  time,
  window,
  log = TRUE,
  se = FALSE,
  ...
)
}
\arguments{
\item{object}{An \code{"\link{egf}"} object.}

\item{what}{A character vector listing one or more variables for which
predicted values are sought.}

\item{time}{A numeric vector supplying time points at which predicted
values are sought. \link{Date} and \link{POSIXt} vectors are
tolerated and coerced to numeric with \code{\link{julian}(time)}.
When \link{time} is missing, time points stored in
\code{\link[=model.frame.egf]{model.frame}(object)} are reused,
and \code{window} is ignored.}

\item{window}{A factor of length \code{length(time)} grouping the elements
of \code{time} by fitting window.
Levels not found in
\code{levels(\link[=model.frame.egf]{model.frame}(object)$window)}
are ignored.}

\item{log}{A logical flag. If \code{FALSE},
then inverse log-transformed predicted values are returned.}

\item{se}{A logical flag. If \code{se = TRUE} and \code{log = TRUE},
then approximate delta method standard errors on predicted values
are reported.
Standard errors are required for subsequent use
of \code{\link{confint.predict.egf}}.
Setting \code{se = TRUE} and \code{log = FALSE} is an error,
as standard errors are not available for inverse log-transformed
predicted values.}

\item{...}{Unused optional arguments.}
}
\value{
A data frame inheriting from class \code{"predict.egf"}, with variables:
\item{var}{
  Predicted variable, from \code{what}.
}
\item{ts}{
  Time series, from
  \code{levels(\link[=model.frame.egf]{model.frame}(object)$ts)}.
}
\item{window}{
  Fitting window, from
  \code{levels(\link[=model.frame.egf]{model.frame}(object)$window)}.
}
\item{time}{
  Time, after possible coercion to numeric.
}
\item{estimate}{
  Predicted value of \code{var} at \code{time} in \code{window},
  conditional on the mixed effects data and fitted model.
}
\item{se}{
  Approximate delta method standard error on \code{estimate}.
  Absent except for calls matching \code{predict(log = TRUE, se = TRUE)}.
}
}
\description{
Computes predicted values of interval incidence, cumulative incidence,
and the per capita growth rate, conditional on observed data and a fitted
nonlinear mixed effects model of epidemic growth.
}
\details{
In the result, \code{estimate[i]} can be interpreted as follows,
for \code{log = FALSE}:
\describe{
\item{\code{interval}}{
  The expected number of cases observed from \code{time[i-1]}
  to \code{time[i]} in \code{window[i]} (interval incidence).
}
\item{\code{cumulative}}{
  The expected number of cases observed up to \code{time[i]}
  in \code{window[i]} (cumulative incidence).
}
\item{\code{rt}}{
  The predicted per capita growth rate at \code{time[i]}.
  This is computed exactly from the differential equation model
  associated with \code{object$model$curve}.
}
}
}
\examples{
object <- egf_cache("egf-1.rds")
zz <- egf_cache("predict-egf-1.rds", predict(object, se = TRUE))
str(zz)

}
