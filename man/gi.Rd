\name{gi}
\alias{gi}
\alias{dgi}
\alias{pgi}
\alias{qgi}
\alias{rgi}
\title{Generation Interval Distribution}
\description{
Generation interval
density function (\code{dgi}), distribution function (\code{pgi}),
quantile function (\code{qgi}), and sampling (\code{rgi}).
Results are conditional on supplied latent and infectious period
distributions.  It is assumed
\itemize{
\item that the latent period and infectious waiting time are independent,
\item that infectiousness is constant over the infectious period, and
\item that the latent and infectious periods are positive and integer-valued
(in arbitrary but common units of time).
}
}
\usage{
dgi(x, latent, infectious)
pgi(q, latent, infectious)
qgi(p, latent, infectious)
rgi(n, latent, infectious)
}
\arguments{
\item{x, q}{
  a numeric vector listing generation intervals.}
\item{p}{
  a numeric vector listing probabilities.}
\item{n}{
  a non-negative integer indicating a sample size.
  If \code{length(n) > 1}, then \code{length(n)} is taken to be the
  sample size.}
\item{latent, infectious}{
  numeric vectors such that \code{latent[i]} and \code{infectious[i]}
  are the probabilities that the latent and infectious periods,
  respectively, are \code{i} units of time.
  It is sufficient to supply probability weights, as internally
  both vectors are divided by their sums.}
}
\value{
A numeric vector with length equal to the that of the first argument
or length \code{n} in the case of \code{rgi}.
}
\references{
Svensson, Å.  (2007).
A note on generation times in epidemic models.
\emph{Mathematical Biosciences},
\emph{208}(1), 300-311.
\doi{10.1016/j.mbs.2006.10.010}
}
\examples{
latent <- c(0.026, 0.104, 0.182, 0.246, 0.318, 0.104,
            0.013, 0.004, 0.003)
m <- length(latent)

infectious <- c(0.138, 0.462, 0.256, 0.078, 0.041, 0.007,
                0.004, 0.004, 0.006, 0.004)
n <- length(infectious)

## Histogram of samples
y <- rgi(1e06, latent, infectious)
hist(y, breaks = seq(0, m + n + 1), freq = FALSE, las = 1,
     ylab = "relative frequency",
     main = "")

## Density and distribution functions
x <- seq(0, m + n + 1, by = 0.02)
fx <- dgi(x, latent, infectious)
Fx <- pgi(x, latent, infectious)
plot(x, fx, type = "l", las = 1, # consistent with histogram
     xlab = "generation interval",
     ylab = "density function")
plot(x, Fx, type = "l", las = 1,
     xlab = "generation interval",
     ylab = "distribution function")

## Quantile function
p <- seq(0, 1, by = 0.001)
qp <- qgi(p, latent, infectious)
plot(p, qp, type = "l", las = 1,
     xlab = "probability",
     ylab = "quantile function")
}
