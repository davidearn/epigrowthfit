\name{timescale}
\alias{timescale}
\title{Compute the Characteristic Time Scale}
\description{
Computes characteristic time scales corresponding to exponential growth
rates.
}
\usage{
timescale(r, units)
}
\arguments{
\item{r}{
  a non-negative numeric vector listing exponential growth rates.}
\item{units}{
  a character string indicating units for the result.
  If missing, then the result is \dQuote{unitless}.}
}
\value{
\code{1/r}, as a \code{\link{difftime}} if \code{units} is not missing.
}
\seealso{
\code{\link{R0}}, \code{\link{finalsize}}.
}
\examples{
r <- 10^(-2:0)
units <- "days"
stopifnot(all.equal(timescale(r), 1 / r),
          all.equal(timescale(r, units), .difftime(1 / r, units)))
}
