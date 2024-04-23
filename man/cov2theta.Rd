\name{cov2theta}
\alias{cov2theta}
\alias{theta2cov}
\title{Compute a Packed Representation of a Covariance Matrix}
\description{
Transform covariances matrices to a \dQuote{packed} representation
or compute the inverse transformation.
}
\usage{
cov2theta(Sigma)
theta2cov(theta)
}
\arguments{
\item{Sigma}{
  an \eqn{n}-by-\eqn{n} real, symmetric positive definite matrix.
  Only the upper triangle is \dQuote{seen}.}
\item{theta}{
  a numeric vector of length \eqn{n(n+1)/2}{n*(n+1)/2} whose first
  \eqn{n} elements are positive.}
}
\details{
An \eqn{n}-by-\eqn{n} real, symmetric, positive definite matrix
\eqn{\Sigma} can be factorized as
\deqn{\Sigma = R' R\,.}{Sigma = R' * R .}
The upper triangular Cholesky factor \eqn{R} can be written as
\deqn{R = R_{1} D^{-1/2} D_{\sigma}^{1/2}\,,}{R = R1 / sqrt(D) * sqrt(Ds) ,}
where
\eqn{R_{1}}{R1} is a unit upper triangular matrix and
\eqn{D = \mathrm{diag}(\mathrm{diag}(R_{1}' R_{1}))}{D = diag(diag(R1' * R1))}
and
\eqn{D_{\sigma} = \mathrm{diag}(\mathrm{diag}(\Sigma))}{Ds = diag(diag(Sigma))}
are diagonal matrices.

\code{cov2theta} takes \eqn{\Sigma} and returns the vector \eqn{\theta}
of length \eqn{n(n+1)/2}{n*(n+1)/2} containing the log diagonal entries
of \eqn{D_{\sigma}} followed by (in column-major order) the strictly
upper triangular entries of \eqn{R_{1}}.  \code{theta2cov} computes the
inverse transformation.
}
