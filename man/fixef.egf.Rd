\name{fixef.egf}
\alias{fixef.egf}
\alias{fixef}
\alias{fixef.egf_no_fit}
\title{Details about Fixed Effect Coefficients}
\description{
Extracts from a model object details about the fixed effect
coefficients, namely segment \code{beta} of the bottom level
parameter vector.
}
\usage{
\method{fixef}{egf}(object, \dots)
}
\arguments{
\item{object}{
  an \code{\link{egf}} object.}
\item{\dots}{
  unused optional arguments.}
}
\value{
A data frame with one row per coefficient and variables:

\item{bottom}{
  label for a bottom level mixed effects model parameter,
  in this case for a fixed effect coefficient.
  This is a string with format \code{"beta[\%d]"}.}
\item{top}{
  name of the top level nonlinear model parameter whose
  fitted value is a function of \code{bottom},
  from \code{\link{egf_top}(object)}.}
\item{term}{
  term from the fixed effects component of the mixed effects
  model formula for parameter \code{top}.}
\item{colname}{
  column name in the fixed effects design matrix
  \code{\link[=model.matrix.egf]{model.matrix}(object, "fixed")}.}
\item{value}{
  coefficient estimate, from segment \code{beta} of
  \code{\link[=coef.egf]{coef}(object, full = TRUE)}.}
}
\seealso{
The generic function \code{\link[nlme]{fixef}}.
}
