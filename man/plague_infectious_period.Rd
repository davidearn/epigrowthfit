\docType{data}
\name{plague_infectious_period}
\alias{plague_infectious_period}
\title{Infectious period distribution (pneumonic plague)}
\format{
A \link[=data.frame]{data frame} with 12 rows and 2 variables:
\describe{
\item{\code{days}}{Infectious period in days, from 1 to 12.}
\item{\code{relfreq}}{Observed relative frequency out of 225 cases.}
}
}
\source{
Obtained from Figure 1B in Gani and Leach (2004).
}
\usage{
data(plague_infectious_period)
}
\description{
Empirical distribution of the infectious period of pneumonic plague.
}
\examples{
data(plague_infectious_period)
plot(relfreq ~ days, data = plague_infectious_period)

}
\references{
Gani R, Leach S. Epidemiological determinants for modeling
pneumonic plague outbreaks. Emerg Infect Dis. 2004;10:608--14.
}
\keyword{datasets}
