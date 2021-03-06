% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prior.R
\name{egf_prior}
\alias{egf_prior}
\alias{Normal}
\title{Prior distributions}
\usage{
Normal(mu = 0, sigma = 1)
}
\arguments{
\item{mu}{Mean.}

\item{sigma}{Standard deviation. Must be positive.}
}
\value{
A \link{list} inheriting from \link{class} \code{"egf_prior"},
with elements:
\item{distribution}{
  A \link{character} string naming the prior distribution.
}
\item{parameters}{
  A named \link{numeric} vector listing parameter values.
}
}
\description{
Functions used by \code{\link{egf}} to specify prior distributions
of nonlinear and dispersion model parameters.
}
