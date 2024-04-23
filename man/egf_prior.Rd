\name{egf_prior}
\alias{egf_prior}
\alias{Normal}
\alias{LKJ}
\alias{Wishart}
\alias{InverseWishart}
\title{Prior Distributions}
\description{
Functions used by \code{\link{egf}} to specify prior distributions
of bottom level mixed effects model parameters.
}
\usage{
Normal(mu = 0, sigma = 1)
LKJ(eta = 1)
Wishart(df, scale, tol = 1e-06)
InverseWishart(df, scale, tol = 1e-06)
}
\arguments{
\item{mu}{
  a numeric vector listing means.}
\item{sigma}{
  a positive numeric vector listing standard deviations.}
\item{eta}{
  a positive numeric vector listing values for the shape parameter,
  with 1 corresponding to a uniform distribution over the space of
  real, symmetric, positive definite matrices with unit diagonal
  elements.  Lesser (greater) values concentrate the probability
  density around such matrices whose determinant is nearer to 0 (1).}
\item{df}{
  a numeric vector listing degrees of freedom.  \code{df} must
  be greater than \code{nrow(scale) - 1}.  If either \code{df} or
  \code{scale} has length greater than 1, then this condition is
  checked elementwise after recycling.}
\item{scale}{
  a list of real, symmetric, positive definite matrices
  or a matrix to be placed in a list of length 1.}
\item{tol}{
  a non-negative number specifying a tolerance for indefiniteness
  of \code{scale}.  All eigenvalues of \code{scale} must exceed
  \code{-tol * rho}, where \code{rho} is the spectral radius of
  \code{scale}.  However, regardless of \code{tol},
  \code{diag(scale)} must be positive, as standard deviations are
  stored on a logarithmic scale.}
}
\value{
A list inheriting from class \code{egf_prior}, with elements:

\item{family}{
  a character string specifying a family of distributions.}
\item{parameters}{
  a named list of numeric vectors specifying parameter values.}
}
\examples{
Normal(mu = 0, sigma = 1)
Normal(mu = -5:5, sigma = c(0.1, 1))

LKJ(eta = 2)

u <- matrix(rnorm(9L), 3L, 3L)
utu <- crossprod(u)
uut <- tcrossprod(u)
Wishart(df = 6, scale = utu)
InverseWishart(df = 6, scale = list(utu, uut))
}
