\name{egf}
\alias{egf}
\alias{egf.egf_model}
\title{Fit Nonlinear Mixed Effects Models of Epidemic Growth}
\description{
Fits nonlinear mixed effects models of epidemic growth
to collections of one or more disease incidence time series.
}
\usage{
egf(model, \dots)

\method{egf}{egf_model}(model,
    formula_ts,
    formula_windows,
    formula_parameters = list(),
    formula_priors = list(),
    data_ts,
    data_windows,
    subset_ts = NULL,
    subset_windows = NULL,
    select_windows = NULL,
    na_action_ts = c("fail", "pass"),
    na_action_windows = c("fail", "omit"),
    control = egf_control(),
    init = list(),
    map = list(),
    fit = TRUE,
    se = FALSE,
    \dots)
}
\arguments{
\item{model}{
  an \R{} object specifying a top level nonlinear model,
  typically of class \code{\link{egf_model}}.}
\item{formula_ts}{
  a formula of the form \code{cbind(time, x) ~ ts}
  specifying one or more disease incidence time series in long format.\cr
  \cr
  \code{ts} must evaluate to a factor (insofar as \code{as.factor(ts)}
  is a factor) grouping the data by time series.
  \code{time} must evaluate to a numeric vector that is increasing
  within levels of \code{ts}.
  \code{\link{Date}}, \code{\link{POSIXct}}, and \code{\link{POSIXlt}}
  vectors are supported and coerced to numeric with
  \code{\link{julian}(time)}.
  Finally, \code{x} must evaluate to a non-negative numeric vector with
  \code{x[i]} equal to the number of cases observed over the interval
  \code{(time[i-1], time[i]]}.
  Edge cases like \code{x[1]} are ignored internally.
  Elements of \code{x} that are not integer-valued are rounded with a
  warning.\cr
  \cr
  \code{formula_ts = cbind(time, x) ~ 1} can be supplied
  when there is only one time series; it is equivalent
  to \code{formula_ts = cbind(time, x) ~ ts} with \code{ts}
  evaluating to \code{rep(factor(1), length(x))}.}
\item{formula_windows}{
  a formula of the form \code{cbind(start, end) ~ ts}
  specifying disjoint fitting windows \code{(start, end]} in long format.
  If \code{formula_ts = cbind(time, x) ~ ts1} and
  \code{formula_windows = cbind(start, end) ~ ts2},
  then observation \code{x[i]} is associated with window
  \code{(start[j], end[j]]} if and only if \code{time[i-1] >= start[j]},
  \code{time[i] <= end[j]}, and \code{ts1[i] == ts2[j]}.}
\item{formula_parameters}{
  a list of formulae of the form \code{parameter ~ terms}
  specifying mixed effects models for top level
  nonlinear model parameters using \pkg{lme4}-like syntax
  (see, e.g., \code{help("lmer", package = "lme4")}).
  Alternatively, a formula of the form \code{~terms} to be recycled
  for all parameters.\cr
  \cr
  A list of parameters for which formulae may be specified can be
  retrieved with \code{\link{egf_top}}.
  Specifically, \code{deparse(parameter)} must be an element of
  \code{\link{egf_top}(model)}.
  The default for parameters not assigned a formula is \code{~1}.}
\item{formula_priors}{
  a list of formulae of the form \code{parameter ~ prior}
  defining priors on:\cr
  (i) top level nonlinear model parameters,\cr
  (ii) fixed effect coefficients and random effect covariance parameters
  (elements of segments \code{beta} and \code{theta} of the bottom level
  parameter vector), or\cr
  (iii) random effect covariance matrices
  (elements of a list \code{Sigma} containing the matrices).\cr
  \cr
  \code{prior} must be a call to a \link[=egf_prior]{prior function}
  with arguments specifying suitable hyperparameters.
  In case (i),
  \code{deparse(parameter)} must be an element of
  \code{\link{egf_top}(model)},
  and hyperparameters supplied on the right hand side must have length 1.
  In cases (ii) and (iii),
  \code{parameter} must be \code{beta}, \code{theta}, or \code{Sigma}
  or a call to \code{[} or \code{[[} referring to a subset or element
  of \code{beta}, \code{theta}, or \code{Sigma}
  (e.g., \code{beta[index]}, where \code{index} is a valid index vector
  for \code{beta}),
  and hyperparameters are recycled to the length of the indicated
  subset.\cr
  \cr
  Expressions \code{prior} and \code{index} are evaluated in the
  corresponding formula environment.}
\item{data_ts, data_windows}{
  data frames, lists, or environments to be searched for variables
  named in the corresponding \code{formula_*} and \code{subset_*}
  arguments.  (\code{formula_parameters} uses \code{data_windows}.)
  Formula environments are searched for variables not found here.}
\item{subset_ts, subset_windows}{
  expressions to be evaluated in the corresponding \code{data_*}
  data frames.
  The value should be a valid index vector for the rows of the data frame.
  Rows that are not indexed are discarded.
  Rows that are indexed are filtered further (e.g., time series
  with zero associated fitting windows are discarded regardless of
  \code{subset_ts}).
  The default is to preserve all rows for further filtering.}
\item{select_windows}{
  an expression indicating additional variables in \code{data_windows}
  (if it is a data frame) to be preserved in the returned object for use
  by methods.
  The default is to preserve nothing.
  A dot \samp{.} is to preserve all variables not occurring in
  \code{formula_windows} or \code{formula_parameters}.
  Outside of these two special cases, evaluation of \code{select}
  follows the implementation of \code{\link{subset.data.frame}}.}
\item{na_action_ts}{
  a character string affecting the handling of \code{NA} in \code{x}
  if \code{formula_ts = cbind(time, x) ~ ts}.
  \code{"fail"} is to throw an error.
  \code{"pass"} is to ignore \code{NA} when fitting and replace
  \code{NA} when predicting.
  \code{NA} in \code{time} and \code{ts} are always an error.}
\item{na_action_windows}{
  a character string affecting the handling of \code{NA} in
  \code{formula_windows} and \code{formula_parameters} variables.
  \code{"fail"} is to throw an error.
  \code{"omit"} is to discard incomplete rows of data.}
\item{control}{
  an \code{\link{egf_control}} object specifying control parameters.}
\item{init}{
  a named list of numeric vectors with possible elements \code{beta},
  \code{theta}, and \code{b}, specifying values to be used in the first
  likelihood evaluation for the so-named segments of the bottom level
  parameter vector.  The default value of each segment is a zero vector,
  with the exception that \code{"(Intercept)"} coefficients in
  \code{beta} have default values computed from supplied time series.
  Use \code{NA} to indicate elements that should retain their default
  value.}
\item{map}{
  a named list of factors with possible elements \code{beta},
  \code{theta}, and \code{b}, each as long as the so-named segment
  of the bottom level parameter.
  Elements of a segment \code{<name>} indexed by
  \code{is.na(map[["<name>"]])}
  are fixed at their initial values, rather than estimated,
  and elements corresponding to a common factor level are constrained
  to have a common value during estimation.
  \code{map[["<name>"]]} can be an index vector for segment
  \code{<name>}, instead of a factor.
  In this case, the indexed elements of that segment are fixed at their
  initial values.}
\item{fit}{
  a logical.  If \code{FALSE}, then \code{egf} returns early
  (\emph{before} fitting) with a partial model object.
  The details of the partial result are subject to change and
  therefore sparsely documented, on purpose \ldots}
\item{se}{
  a logical.  If \code{TRUE}, then the Hessian matrix of the negative
  log marginal likelihood function is computed and inverted to
  approximate the joint covariance matrix of segments \code{beta} and
  \code{theta} of the bottom level parameter vector.
  Standard errors on the fitted values of all top level nonlinear model
  parameters are computed approximately using the delta method.
  Computations are preserved in the model object for reuse by methods.}
\item{\dots}{
  additional arguments passed from or to other methods.}
}
\value{
A list inheriting from class \code{egf}.
See topic \code{\link{egf-class}} for class documentation.
}
\details{
Users attempting to set arguments \code{formula_priors}, \code{init},
and \code{map} should know the structure of the bottom level parameter
vector.  It is described under topic \code{\link{egf-class}}.

If
\preformatted{formula_ts = cbind(time, x) ~ ts1
formula_windows = cbind(start, end) ~ ts2}
then it is expected that \code{time}, \code{start}, and \code{end}
(after coercion to numeric) measure time on the same scale.
To be precise, numeric times should have a common unit of measure
and, at least within time series, represent displacements from a
common reference time.
These conditions will always hold if \code{time}, \code{start}, and
\code{end} all evaluate to \code{\link{Date}}, \code{\link{POSIXct}},
or \code{\link{POSIXlt}} vectors.

When day of week effects are estimated, numeric times are interpreted
as numbers of days since midnight on January 1, 1970, so that time
points can be mapped unambiguously to days of week.
Furthermore, in this case, \code{time} (after coercion to numeric) is
required to be integer-valued with one day spacing in all time series.
This means that
\preformatted{isTRUE(all.equal(time, round(time))) && all(diff(round(time)) == 1)}
must be \code{TRUE} in each level of \code{ts1}.
These conditions ensure that intervals between successive time points
represent exactly one day of week.
}
\seealso{
The many methods for class \code{egf},
listed by \code{\link{methods}(class = "egf")}.

Function \code{\link{egf1}}, an interface to \code{egf} providing
simpler usage for the case in which there is exactly one time series and
exactly one fitting window.
}
\examples{
## Simulate 'N' incidence time series exhibiting exponential growth
set.seed(180149L)
N <- 10L
f <- function(time, r, c0) {
    lambda <- diff(exp(log(c0) + r * time))
    c(NA, rpois(lambda, lambda))
}
time <- seq.int(0, 40, 1)
r <- rlnorm(N, -3.2, 0.2)
c0 <- rlnorm(N, 6, 0.2)
data_ts <-
    data.frame(country = gl(N, length(time), labels = LETTERS[1:N]),
               time = rep.int(time, N),
               x = unlist(Map(f, time = list(time), r = r, c0 = c0)))
rm(f, time)

## Define fitting windows (here, two per time series)
data_windows <-
    data.frame(country = gl(N, 1L, 2L * N, labels = LETTERS[1:N]),
               wave = gl(2L, 10L),
               start = c(sample(seq.int(0, 5, 1), N, TRUE),
                         sample(seq.int(20, 25, 1), N, TRUE)),
               end = c(sample(seq.int(15, 20, 1), N, TRUE),
                       sample(seq.int(35, 40, 1), N, TRUE)))

## Estimate the generative model
m1 <-
    egf(model = egf_model(curve = "exponential", family = "pois"),
        formula_ts = cbind(time, x) ~ country,
        formula_windows = cbind(start, end) ~ country,
        formula_parameters = ~(1 | country:wave),
        data_ts = data_ts,
        data_windows = data_windows,
        se = TRUE)

## Re-estimate the generative model with:
## * Gaussian prior on beta[1L]
## * LKJ prior on all random effect covariance matrices
##   (here there happens to be just one)
## * initial value of 'theta' set explicitly
## * theta[3L] fixed at initial value
m2 <-
    update(m1,
           formula_priors = list(beta[1L] ~ Normal(mu = -3, sigma = 1),
                                 Sigma ~ LKJ(eta = 2)),
           init = list(theta = c(log(0.5), log(0.5), 0)),
           map = list(theta = 3L))
}
