\name{egf}
\alias{egf}
\alias{egf.egf_model}
\title{Fit nonlinear mixed effects models of epidemic growth}
\usage{
egf(model, ...)

\method{egf}{egf_model}(
  model,
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
  fit = TRUE,
  se = FALSE,
  init = list(),
  map = list(),
  ...
)
}
\arguments{
\item{model}{An \code{"\link{egf_model}"} object specifying a top level
nonlinear model to be estimated.}

\item{...}{Arguments passed to methods by the generic function.}

\item{formula_ts}{A \link{formula} of the form \code{cbind(time, x) ~ ts}
specifying one or more disease incidence time series in long format.
\code{ts} must evaluate to a \link{factor}
(insofar as \code{\link{as.factor}(ts)} is a factor)
grouping the data by time series.
\code{time} must evaluate to a \link{numeric} vector
that is increasing within \link{levels} of \code{ts}.
\link{Date} and \link{POSIXt} vectors are tolerated and
coerced to numeric with \code{\link{julian}(time)}.
Finally, \code{x} must evaluate to a non-negative numeric vector
with \code{x[i]} equal to the number of cases observed over the
interval \code{(time[i-1], time[i]]}.
Edge cases like \code{x[1]} are ignored internally.
Nonintegral elements of \code{x} are rounded with a warning.
\code{formula_ts = cbind(time, x) ~ 1} can be supplied
when there is only one time series; it is equivalent
to \code{formula_ts = cbind(time, x) ~ ts} with \code{ts}
evaluating to \code{\link{rep}(factor(1), length(x))}.}

\item{formula_windows}{A \link{formula} of the form \code{cbind(start, end) ~ ts}
specifying disjoint fitting windows \code{(start, end]} in long format.
If
\code{formula_ts = cbind(time, x) ~ ts1}
and
\code{formula_windows = cbind(start, end) ~ ts2},
then observation \code{x[i]}
is associated with window \code{(start[j], end[j]]} if and only if
\code{time[i-1] >= start[j]},
\code{time[i] <= end[j]}, and
\code{ts1[i] == ts2[j]}.}

\item{formula_parameters}{A \link{list} of \link{formula}e of the form \code{parameter ~ terms}
specifying mixed effects models for top level nonlinear model parameters,
using \pkg{lme4}-like syntax (see \code{?lme4::lmer}).
Alternatively, a formula of the form \code{~terms} to be recycled for
all parameters.
A list of parameters for which formulae may be specified can be retrieved
with \code{\link{egf_top}}.
Specifically, \code{\link{deparse}(parameter)} must be an element of
\code{\link{egf_top}(model)}.
The default for parameters not assigned a formula is \code{~1}.}

\item{formula_priors}{A \link{list} of \link{formula}e of the form \code{parameter ~ prior}
defining priors on:\cr
(i) top level nonlinear model parameters,\cr
(ii) fixed effects coefficients and random effect covariance parameters
(elements of segments \code{beta} and \code{theta} of the bottom level
parameter vector, or\cr
(iii) random effect covariance matrices
(elements of a \link{list} \code{Sigma} containing the matrices).\cr
\code{prior} must be a \link{call} to a \link[=egf_prior]{prior} function
with arguments specifying suitable hyperparameters.
In case (i),
\code{deparse(parameter)} must be an element of
\code{\link{egf_top}(model)},
and hyperparameters supplied on the right hand side must have length 1.
In cases (ii) and (iii),
\code{parameter} must be \code{beta}, \code{theta}, or \code{Sigma}
or a call to \code{\link{[}} or \code{\link{[[}} referring to a subset
or element of \code{beta}, \code{theta}, or \code{Sigma}
(e.g., \code{beta[index]}, where \code{index} is a valid index vector
for \code{beta}),
and hyperparameters are recycled to the length of the indicated subset.
All expressions \code{prior} and \code{index} are evaluated in the
corresponding formula environment.}

\item{data_ts, data_windows}{\link[=data.frame]{Data frame}s, \link{list}s, or \link{environment}s
to be searched for variables named in the corresponding formulae and
subset expressions. (\code{data_windows} is also searched for variables
named in \code{formula_parameters}.) Formula environments are searched
for variables not found here.}

\item{subset_ts, subset_windows}{Expressions to be evaluated in the corresponding data frame.
The result should be a valid index vector for the rows of the data frame
(see \code{\link{[.data.frame}}).
Rows that are not indexed are discarded.
Rows that are indexed are filtered further
(e.g., time series with zero associated fitting windows are discarded
regardless of \code{subset_ts}).
The default (\code{\link{NULL}}) is to preserve all rows for further
filtering.}

\item{select_windows}{An expression indicating additional variables
in \code{data_windows} to be preserved in the returned object for use
by methods.
Usage requires that \code{data_windows} is a \link[=data.frame]{data frame}.
The default (\code{\link{NULL}}) is to preserve nothing.
A dot \samp{.} is to preserve all variables not occurring in
\code{formula_windows} or \code{formula_parameters}.
Outside of these two special cases, the expression is evaluated
similarly to argument \code{select} of function \code{\link{subset}}.}

\item{na_action_ts}{A \link{character} string affecting the handling of \code{\link{NA}}
in \code{x} if \code{formula_ts = cbind(time, x) ~ ts}.
\code{"fail"} is to throw an error.
\code{"pass"} is to ignore \code{NA} when fitting and replace \code{NA}
when predicting.
Note that \code{NA} in \code{time} and \code{ts} are always an error.}

\item{na_action_windows}{A \link{character} string affecting the handling of \code{\link{NA}}
in \code{formula_windows} and \code{formula_parameters} variables.
\code{"fail"} is to throw an error.
\code{"omit"} is to discard incomplete rows of data.}

\item{control}{An \code{"\link{egf_control}"} object specifying control parameters.}

\item{fit}{A \link{logical} flag. If \code{FALSE}, then \code{egf} returns early
(\emph{before} fitting) with a partial model object.}

\item{se}{A \link{logical} flag.
If \code{TRUE}, then the Hessian matrix of the negative log likelihood
function is computed and inverted to approximate the joint covariance
matrix of segments \code{beta} and \code{theta} of the bottom level
parameter vector.
In addition, the standard errors of the fitted values of all top level
nonlinear model parameters are computed approximately using the delta
method.
Computations are preserved in the model object for reuse by methods.}

\item{init}{A named \link{list} of \link{numeric} vectors with possible elements
\code{beta}, \code{theta}, and \code{b}, specifying values to be used
in the first likelihood evaluation for the so-named segments of the
bottom level parameter vector.
The default value of each segment is a zero vector, with the exception
that \code{"(Intercept)"} coefficients in \code{beta} have default values
computed from supplied time series.
To specify only a subset of a segment, use \code{\link{NA}} to indicate
elements that should retain their default value.}

\item{map}{A named list of \link{factor}s with possible elements \code{beta},
\code{theta}, and \code{b}, each as long as the so-named segment
of the bottom level parameter.
Elements of segment \code{name} indexed by \code{is.na(map[["name"]])},
are fixed at their initial values, rather than estimated, and elements
corresponding to a common factor level are are constrained to have a
common value during estimation.
\code{map[["name"]]} can be an index vector for segment \code{name},
instead of a factor. In this case, the elements of segment \code{name}
indexed by \code{map[["name"]]} are fixed at their initial values.}
}
\value{
A \link{list} inheriting from \link{class} \code{"egf"}
or \code{"egf_no_fit"} depending on the value of \code{fit}.
See topic \code{\link{egf-class}} for class documentation.
}
\description{
Fits nonlinear mixed effects models of epidemic growth
to collections of one or more disease incidence time series.
}
\details{
Users attempting to set arguments \code{formula_priors}, \code{init}, and
\code{map} should know the structure of the bottom level parameter vector.
It is described under topic \code{\link{egf-class}}.

If
\code{formula_ts = cbind(time, x) ~ ts1}
and
\code{formula_windows = cbind(start, end) ~ ts2},
then
it is expected that \code{time}, \code{start}, and \code{end}
(after possible coercion to \link{numeric}) measure time
on the same scale. To be precise, numeric times should have
a common unit of measure and, at least within time series,
represent displacements from a common reference time.
These conditions will always hold if \code{time}, \code{start},
and \code{end} all evaluate to \link{Date} or \link{POSIXt}
vectors.

When day of week effects are estimated (model$day_of_week > 0),
numeric times must be interpretable as numbers of days since
\code{1970-01-01 00:00:00}, so that time points can be mapped
unambiguously to days of week. Furthermore, in this case,
\code{time} (after possible coercion to \link{numeric})
is required to be integer-valued with one day spacing in all
time series. This means that
\code{\link{all.equal}(time, \link{round}(time))}
and
\code{\link{all}(\link{diff}(\link{round}(time)) == 1)}
must both be \code{TRUE} in each level of \code{ts1}.
These conditions ensure that intervals between successive
time points each represent exactly one day of week.
}
\examples{
## Simulate 'N' incidence time series exhibiting exponential growth
set.seed(180149L)
N <- 10L
time <- seq.int(0, 40, 1)
mu <- c(-3.2, 6)
sigma <- c(0.2, 0.2)
r <- rlnorm(N, -3.2, 0.2)
c0 <- rlnorm(N, 6, 0.2)
f <- function(time, r, c0) {
    lambda <- diff(exp(log(c0) + r * time))
    c(NA, rpois(lambda, lambda))
}
data_ts <-
    data.frame(country = gl(N, length(time), labels = LETTERS[1:N]),
               time = rep.int(time, N),
               x = unlist(Map(f, time = list(time), r = r, c0 = c0)))

## Define fitting windows (here, two per time series)
data_windows <-
    data.frame(country = gl(N, 1L, 2L * N, labels = LETTERS[1:N]),
               wave = gl(2L, 10L),
               start = c(sample(seq.int(0, 5, 1), N, TRUE),
                         sample(seq.int(20, 25, 1), N, TRUE)),
               end = c(sample(seq.int(15, 20, 1), N, TRUE),
                       sample(seq.int(35, 40, 1), N, TRUE)))

## Estimate the generative model
m1 <- egf_cache("egf-1.rds", {
    egf(model = egf_model(curve = "exponential", family = "pois"),
        formula_ts = cbind(time, x) ~ country,
        formula_windows = cbind(start, end) ~ country,
        formula_parameters = ~(1 | country:wave),
        data_ts = data_ts,
        data_windows = data_windows,
        se = TRUE)
})

## Re-estimate the generative model with:
## * Gaussian prior on 'beta[1L]'
## * LKJ prior on all random effect covariance matrices (here there happens to be just one)
## * initial value of 'theta' set explicitly
## * 'theta[3L]' fixed at initial value
m2 <- egf_cache("egf-2.rds", {
    update(m1,
           formula_priors = list(beta[1L] ~ Normal(mu = -3, sigma = 1),
                                 Sigma ~ LKJ(eta = 2)),
           init = list(theta = c(log(0.5), log(0.5), 0)),
           map = list(theta = 3L))
})

}
