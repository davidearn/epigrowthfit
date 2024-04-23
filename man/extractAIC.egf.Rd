\name{extractAIC.egf}
\alias{extractAIC.egf}
\title{Extract the (Generalized) AIC}
\description{
Extracts from a model object the generalized Akaike Information Criterion
(AIC).
}
\usage{
\method{extractAIC}{egf}(fit, scale, k = 2, \dots)
}
\arguments{
\item{fit}{
  an \code{\link{egf}} object.}
\item{scale}{
  unused argument, for generic consistency.}
\item{k}{
  a number giving a weight for the equivalent degrees of freedom.
  \code{k=2} and \code{k=log(nobs(fit))} give the standard
  Akaike Information Criterion and Bayesian Information Criterion.}
\item{\dots}{
  unused optional arguments.}
}
\value{
An numeric vector of length 2 giving the equivalent degrees of freedom
and criterion value.
}
\seealso{
The generic function \code{\link{extractAIC}}.
}
