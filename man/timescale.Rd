\name{timescale}
\alias{timescale}
\title{Compute the characteristic time scale}
\usage{
timescale(r, units)
}
\arguments{
\item{r}{A non-negative numeric vector listing exponential growth rates.}

\item{units}{A character string indicating units for the result.
If missing, then the result is unitless.}
}
\value{
\code{1 / r}, as a \code{\link{difftime}} if \code{units} is not missing.
}
\description{
Computes characteristic time scales corresponding to exponential growth rates.
}
\examples{
r <- 10^(-2:0)
units <- "days"
stopifnot(exprs = {
	all.equal(timescale(r), 1 / r)
	all.equal(timescale(r, units), .difftime(1 / r, units))
})

}
\seealso{
\code{\link{R0}}, \code{\link{finalsize}}
}
