\name{timescale}
\alias{timescale}
\title{Compute doubling time}
\usage{
timescale(r, per = NULL)
}
\arguments{
\item{r}{A non-negative numeric vector listing exponential growth rates.}

\item{per}{A positive number indicating that \code{r} is a rate per
\code{per} days, in which case the result is printed with units.
Use the default (\code{NULL}) if \code{r} is unitless.}
}
\value{
\code{log(2) / r} with \code{"tdoubling"} prepended to its \code{class}
attribute.
\code{per} is retained as an attribute for use by \code{print.tdoubling}.
}
\description{
Computes doubling times corresponding to exponential growth rates.
}
\examples{
r <- 10^(-2:0)
tdoubling <- timescale(r)
all.equal(tdoubling, log(2) / r)

## 'per' affects printing
for (i in c(1, 2, 7, 8, 365, 366)) {
    attr(tdoubling, "per") <- i
    print(tdoubling)
}

}
\seealso{
Other epidemic parameters: 
\code{\link{R0}()},
\code{\link{finalsize}()}
}
\concept{epidemic parameters}
