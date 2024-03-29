\name{egf_prior}
\alias{egf_prior}
\alias{Normal}
\alias{LKJ}
\alias{Wishart}
\alias{InverseWishart}
\title{Prior distributions}
\usage{
Normal(mu = 0, sigma = 1)

LKJ(eta = 1)

Wishart(df, scale, tol = 1e-06)

InverseWishart(df, scale, tol = 1e-06)
}
\arguments{
\item{mu}{A numeric vector listing means.}

\item{sigma}{A positive numeric vector listing standard deviations.}

\item{eta}{A positive numeric vector listing values for the shape parameter,
with 1 corresponding to a uniform distribution over the space of
symmetric positive definite matrices with unit diagonal elements.
Lesser (greater) values concentrate the probability density around
such matrices whose determinant is nearer to 0 (1).}

\item{df}{A numeric vector listing degrees of freedom.
\code{df} must be greater than \code{nrow(scale) - 1}.
(If either \code{df} or \code{scale} has length greater than 1,
then this condition is checked pairwise after recycling.)}

\item{scale}{A list of symmetric positive definite \link{numeric} matrices,
or a matrix to be placed in a list of length 1.}

\item{tol}{A non-negative number specifying a tolerance for non-positive definiteness
of \code{scale}. All eigenvalues of \code{scale} must exceed
\code{-tol * rho}, where \code{rho} is the spectral radius of \code{scale}.
(However, regardless of \code{tol}, \code{\link{diag}(scale)} must
be positive, as standard deviations are processed on the log scale.)}
}
\value{
A list inheriting from class \code{"egf_prior"}, with elements:
\item{family}{
  A character string naming a family of distributions.
}
\item{parameters}{
  A named list of numeric vectors specifying parameter values.
}
}
\description{
Functions used by \code{\link{egf}} to specify prior distributions
of bottom level mixed effects model parameters.
}
\examples{
Normal(mu = 0, sigma = 1)
Normal(mu = -5:5, sigma = c(0.1, 1))

LKJ(eta = 2)

U <- matrix(rnorm(9L), 3L, 3L)
UTU <- t(U) \%*\% U
UUT <- U \%*\% t(U)
Wishart(df = 6, scale = UTU)
InverseWishart(df = 6, scale = list(UTU, UUT))

}
