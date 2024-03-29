\name{R0}
\alias{R0}
\title{Compute the basic reproduction number}
\usage{
R0(r, breaks, probs)
}
\arguments{
\item{r}{A non-negative numeric vector listing initial exponential growth rates.}

\item{breaks}{An increasing numeric vector of length 2 or greater listing
break points in the support of the generation interval distribution,
in reciprocal units of \code{r}.}

\item{probs}{A numeric vector of length \code{length(breaks)-1}.
\code{probs[i]} is the probability that the generation interval
is between \code{breaks[i]} and \code{breaks[i+1]}.
If \code{sum(probs) != 1}, then \code{probs} is replaced
with \code{probs / sum(probs)}.}
}
\value{
A numeric vector of length \code{length(r)} listing
basic reproduction numbers.
}
\description{
Computes the basic reproduction number as a function of the
initial exponential growth rate, conditional on a binned
generation interval distribution.
}
\section{Computation}{

For an initial exponential growth rate \code{r},
the basic reproduction number \code{R0} is computed as

\code{R0(r) = r / sum(probs * (exp(-r * breaks[-n]) - exp(-r * breaks[-1L])) / (breaks[-1L] - breaks[-n]))},

where \code{n = length(breaks)}.
}

\examples{
r <- seq(0, 1, 0.02)
breaks <- seq(0, 20, 1)
probs <- diff(pgamma(breaks, shape = 1, scale = 2.5))

plot(r, R0(r, breaks, probs), las = 1,
     xlab = "initial exponential growth rate",
     ylab = "basic reproduction number")

}
\references{
Wallinga J, Lipsitch M. How generation intervals shape the relationship
between growth rates and reproductive numbers. Proc R Soc Lond B Biol Sci.
2007;274:599--604.
}
\seealso{
\code{\link{timescale}}, \code{\link{finalsize}}
}
