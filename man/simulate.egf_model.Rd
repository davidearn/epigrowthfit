\name{simulate.egf_model}
\alias{simulate.egf_model}
\alias{coef.simulate.egf_model}
\alias{getCall.simulate.egf_model}
\title{Simulate Incidence Time Series}
\description{
Simulates equally spaced incidence time series according to a specified
nonlinear model.  Top level nonlinear model parameters vary between
time series according to a fixed intercept model \code{~ts} or random
intercept model \code{~(1 | ts)}.
}
\usage{
\method{simulate}{egf_model}(object, nsim = 1, seed = NULL,
         mu, Sigma = NULL, tol = 1e-06, cstart = 0, tend = 100, \dots)
}
\arguments{
\item{object}{an \code{\link{egf_model}} object specifying a top level
  nonlinear model to be simulated.}
\item{nsim}{a positive integer indicating a number of time series.}
\item{seed}{an integer used to set the \link{RNG} state before
  simulation or, otherwise, \code{NULL}; see \code{\link{simulate}}.}
\item{mu}{a numeric vector listing means across time series of top
  level nonlinear model parameters (link scale).  It is assumed that
  elements are ordered as \code{\link{egf_top}(object)}.}
\item{Sigma}{a real, symmetric positive definite matrix to be used
  as the covariance matrix corresponding to \code{mu}.  The default
  is equivalent to a zero matrix and is handled specially.}
\item{tol}{a non-negative number indicating a tolerance for indefinite
  \code{Sigma}.  Eigenvalues of \code{Sigma} must exceed \code{-tol}
  times its spectral radius.  \code{diag(Sigma)} must be positive
  regardless of \code{tol}, as standard deviations are handled on a
  log scale.}
\item{cstart}{a number indicating a threshold value of cumulative
  incidence.  Left endpoints of suggested fitting windows are those
  times when cumulative incidence first exceeds this threshold.}
\item{tend}{a positive number.  Simulated time series run from time 0
  to time \code{tend} with unit spacing.  For nonlinear models of
  expected cumulative incidence with inflections, this argument is
  ignored and set to \code{tinfl+1}, where \code{tinfl} is the time
  of inflection.}
\item{\dots}{unused optional arguments.}
}
\value{
A list inheriting from class \code{simulate.egf_model}, with elements:

\item{model}{copy of argument \code{object}.}
\item{formula_ts}{a formula, always \code{cbind(time, x) ~ ts},
  expressing how simulated time series are stored in \code{data_ts}.}
\item{formula_windows}{a formula, always \code{cbind(start, end) ~ ts},
  expressing how fitting window endpoints are stored in \code{data_windows}.}
\item{formula_parameters}{a formula specifying the generative model.
  If \code{Sigma = NULL}, then the formula is \code{~1}
  if \code{nsim = 1} and \code{~ts} if \code{nsim > 1}.
  Otherwise, it is \code{~(1 | ts)}.}
\item{data_ts}{a data frame with variables \code{ts}, \code{time}, and
  \code{x} storing \code{nsim} simulated time series in long format.}
\item{data_windows}{a data frame with \code{nsim} rows and variables
  \code{ts}, \code{start}, and \code{end} suggesting a fitting window
  for each simulated time series.  Start times are determined by
  \code{cstart}.  End times are always the last time point in the
  corresponding time series.}
\item{init}{a named list of the form \code{list(beta, theta)}
  giving the full bottom level parameter vector of the generative
  model.}
\item{Y}{a numeric matrix with \code{nsim} rows and \code{length(mu)}
  columns listing top level nonlinear model parameter values for each
  time series.  If \code{Sigma = NULL}, then the row vectors of are
  all \code{mu}.  Otherwise, \code{Y} is (conceptually) the result of
  \code{MASS::mvrnorm(nsim, mu, Sigma, tol)}.}
\item{call}{the call to \code{\link{simulate}}, allowing for updates
  to the \code{simulate.egf_model} object via \code{\link{update}}.}

Attribute \code{RNGstate} preserves the RNG state prior to simulation,
making the result reproducible.
}
\seealso{
\code{\link{egf.simulate.egf_model}}.
}
\examples{
r <- 0.04
c0 <- 400
s <- 0.2

mu <- log(c(r, c0))
Sigma <- diag(rep.int(s^2, length(mu)))

zz <- egf_cache("simulate-egf_model-1.rds",
    simulate(object = egf_model(curve = "exponential", family = "pois"),
             nsim = 20L,
             seed = 202737L,
             mu = mu,
             Sigma = Sigma,
             cstart = 10))
str(zz)
}
