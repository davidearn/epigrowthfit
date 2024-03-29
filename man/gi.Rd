\name{gi}
\alias{gi}
\alias{dgi}
\alias{pgi}
\alias{qgi}
\alias{rgi}
\title{Generation interval distribution}
\usage{
dgi(x, latent, infectious)

pgi(q, latent, infectious)

qgi(p, latent, infectious)

rgi(n, latent, infectious)
}
\arguments{
\item{x, q}{A numeric vector listing generation intervals.}

\item{latent, infectious}{Numeric vectors such that \code{latent[i]} and \code{infectious[i]}
are the probabilities that the latent and infectious periods,
respectively, are \code{i} units of time.
It is sufficient to supply probability weights, as both vectors
are divided by their sums internally.}

\item{p}{A numeric vector listing probabilities.}

\item{n}{A non-negative integer indicating a sample size.
If \code{length(n) > 1}, then \code{length(n)} is taken
to be the sample size.}
}
\value{
A numeric vector with length equal to the that of the first argument,
or length \code{n} in the case of \code{rgi}.
}
\description{
Generation interval
density function (\code{dgi}), distribution function (\code{pgi}),
quantile function (\code{qgi}), and sampling (\code{rgi}).
Results are conditional on supplied latent and infectious period
distributions. It is assumed
\itemize{
\item that the latent period and infectious waiting time are independent,
\item that infectiousness is constant over the infectious period, and
\item that the latent and infectious periods are positive and integer-valued
(in arbitrary but like units of time).
}
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
\references{
Svensson, Å. A note on generation times in epidemic models.
Math Biosci. 2007;208:300--11.
}
