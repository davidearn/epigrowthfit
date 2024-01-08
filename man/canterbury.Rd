\docType{data}
\name{canterbury}
\alias{canterbury}
\title{Canterbury wills}
\format{
A \link[=data.frame]{data frame} with 32681 rows and 4 variables:
\describe{
\item{\code{Date}}{
  A \link{Date} vector. Spacing varies, as only nonzero counts of wills
  are included.
}
\item{\code{wills}}{
  An \link{integer} vector. \code{wills[i]} is the number of wills written
  on \code{Date[i]}.
}
\item{\code{outbreak}}{
  An ordered \link{factor}, \code{split}ting the time series by plague
  outbreak. \link[=levels]{Levels} indicate the years in which outbreaks
  took place: 1563, 1578, ..., 1647, and 1665 (the Great Plague of London).
}
\item{\code{severity}}{
  An ordered \link{factor}, classifying outbreaks as \code{"minor"}
  or \code{"major"}. An outbreak is classified as major if and only if
  plague deaths per week per 1000 individuals exceeded 5 at least once.
}
}
}
\source{
These data were retrieved from the National Archives (UK) in 2018
using this
\href{https://www.nationalarchives.gov.uk/help-with-your-research/research-guides/wills-1384-1858/}{research guide}.
}
\usage{
data(canterbury)
}
\description{
Counts of wills probated in the Prerogative Court of Canterbury
during 24 plague epidemics in 16th and 17th century London.
}
\examples{
data(canterbury)
subset(canterbury, outbreak == 1665, -outbreak) # Great Plague of London
}
\keyword{datasets}
