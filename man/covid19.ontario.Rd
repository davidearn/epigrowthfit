\name{covid19.ontario}
\alias{covid19.ontario}
\docType{data}
\keyword{datasets}
\title{COVID-19 in Ontario, Canada}
\description{
Time series of COVID-19 cases and tests in Ontario, Canada,
daily from February 8, 2020 to May 1, 2022.
}
\usage{
data(covid19.ontario, package = "epigrowthfit")
}
\format{
A \link[=data.frame]{data frame} with 814 rows and 3 variables:
\describe{
\item{\code{date}}{a \link{Date} vector.}
\item{\code{cases}}{an \link{integer} vector.  \code{cases[i]} is
  the number of cases confirmed by Ontario public health units prior
  to \code{date[i]}.  This number includes resolved and fatal cases
  as well as reinfections.}
\item{\code{tests}}{an \link{integer} vector.  \code{tests[i]} is
  the number of tests completed prior to \code{date[i]}.  This number
  includes repeated tests by individuals except prior to April 15,
  2020, when individuals were counted at most once.}
}
}
\source{
This data set is a processed subset of a larger data set downloaded
on 2024-01-10 from the \href{https://data.ontario.ca/}{Ontario Data Catalogue}.
It is updated using an installed script; see
\code{\link{system.file}("scripts", "covid19.ontario.R", package = "epigrowthfit")}.
}
\examples{
data(covid19.ontario, package = "epigrowthfit")
plot(1 + diff(c(NA, cases)) ~ date, data = covid19.ontario, log = "y")
}
