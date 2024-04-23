\name{R0}
\alias{R0}
\title{Compute the Basic Reproduction Number}
\description{
Computes the basic reproduction number as a function of the initial
exponential growth rate, conditional on a binned generation interval
distribution.
}
\usage{
R0(r, breaks, probs)
}
\arguments{
\item{r}{
  a non-negative numeric vector listing initial exponential growth
  rates.}
\item{breaks}{
  an increasing numeric vector of length 2 or greater listing break
  points in the support of the generation interval distribution,
  in reciprocal units of \code{r}.}
\item{probs}{
  a numeric vector of length \code{length(breaks)-1}.
  \code{probs[i]} is the probability that the generation interval
  is between \code{breaks[i]} and \code{breaks[i+1]}.
  It is sufficient to supply probability weights, as internally
  the vector is divided by its sum.}
}
\value{
A numeric vector listing basic reproduction numbers.
}
\section{Computation}{
For an initial exponential growth rate \code{r},
the basic reproduction number is computed as
\preformatted{r / sum(probs * (exp(-r * breaks[-n]) - exp(-r * breaks[-1L])) / (breaks[-1L] - breaks[-n])) ,}
where \code{n = length(breaks)}.
}
\seealso{
\code{\link{timescale}}, \code{\link{finalsize}}.
}
\references{
Wallinga, J. & Lipsitch M. (2007).
How generation intervals shape the relationship between growth rates
and reproductive numbers.
\emph{Proceedings of the Royal Society B: Biological Sciences},
\emph{274}(1609), 599-604.
\doi{10.1098/rspb.2006.3754}
}
\examples{
r <- seq(0, 1, 0.02)
breaks <- seq(0, 20, 1)
probs <- diff(pgamma(breaks, shape = 1, scale = 2.5))

plot(r, R0(r, breaks, probs), las = 1,
     xlab = "initial exponential growth rate",
     ylab = "basic reproduction number")
}
