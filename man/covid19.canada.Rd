\docType{data}
\name{covid19.canada}
\alias{covid19.canada}
\title{COVID-19 in Canadian provinces and territories}
\format{
A \link[=data.frame]{data frame} with 9881 rows and 3 variables:
\describe{
\item{\code{province}}{
  A \link{factor}. Canadian province or territory (postal abbreviation).
}
\item{\code{date}}{
  A \link{Date} vector. Daily within each level of \code{province},
  except prior to invasion, when less frequent reports were common.
}
\item{\code{confirmed}}{
  An \link{integer} vector. \code{confirmed[i]} is the number of cases
  confirmed up to the start of \code{date[i]} in \code{province[i]}.
}
}
}
\source{
This data set is a processed subset of a larger data set downloaded
on 2024-01-10 from Michael Li's public
\href{https://github.com/wzmli/COVID19-Canada}{GitHub repository}.
It is updated using an installed script; see
\code{\link{system.file}("scripts", "covid19.canada.R", package = "epigrowthfit")}.
}
\usage{
data(covid19.canada, package = "epigrowthfit")
}
\description{
Daily confirmations of COVID-19 in Canadian provinces
and territories, from first confirmation to May 1, 2022.
}
\examples{
data(covid19.canada, package = "epigrowthfit")
subset(covid19.canada, province == "ON", -province)

}
\keyword{datasets}
