% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{finalsize}
\alias{finalsize}
\title{expected final size of an epidemic}
\usage{
finalsize(R0)
}
\arguments{
\item{R0}{basic reproduction number}
}
\description{
For a given R0, return the proportion of the population
expected to be infected if an epidemic is seeded and
all individuals are initially susceptible.  The standard
final size formula (Kermack and McKendrick 1927,
Ma and Earn 2006) is used.
}
\references{
\insertRef{KermMcKe27}{epigrowthfitPNAS}

\insertRef{MaEarn06}{epigrowthfitPNAS}
}