\docType{data}
\name{covid19.canada}
\alias{covid19.canada}
\title{COVID-19 in Canadian provinces and territories}
\format{
A \link[=data.frame]{data frame} with 5451 rows and 4 variables:
\describe{
\item{\code{province}}{
  A \link{factor}. Canadian province or territory (postal abbreviation).
}
\item{\code{Date}}{
  A \link{Date} vector. Daily within each level of \code{province},
  except prior to invasion, when less frequent reports were common.
}
\item{\code{cases_tot}}{
  An \link{integer} vector. \code{cases_tot[i]} is the number of cases
  confirmed up to the end of \code{Date[i]} in \code{province[i]}.
}
\item{\code{cases_new}}{
  An \link{integer} vector. \code{cases_new[i]} is the number of cases
  confirmed from the end of \code{Date[i-1]} to the end of \code{Date[i]}
  in \code{province[i]}.
}
}
}
\source{
Raw data were downloaded from Michael Li's public
\href{https://github.com/wzmli}{Github repository}.
Up-to-date data can be downloaded
\href{https://wzmli.github.io/COVID19-Canada/COVID19_Canada.csv}{here}.
}
\usage{
data(covid19.canada)
}
\description{
Daily confirmations of COVID-19 in Canadian provinces
and territories, from first confirmation to May 8, 2021.
}
\examples{
data(covid19.canada)
subset(covid19.canada, province == "ON", -province)

}
\keyword{datasets}
