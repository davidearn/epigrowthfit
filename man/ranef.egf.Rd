\name{ranef.egf}
\alias{ranef.egf}
\alias{ranef}
\alias{ranef.egf_no_fit}
\title{Details about Random Effect Coefficients}
\description{
Extracts from a model object details about the random effect
coefficients, namely segment \code{b} of the bottom level
parameter vector.
}
\usage{
\method{ranef}{egf}(object, makeSigma = FALSE, \dots)
}
\arguments{
\item{object}{
  an \code{\link{egf}} object.}
\item{makeSigma}{
  a logical flag.  If \code{TRUE}, then random effect covariance matrices
  are constructed from segment \code{theta} of
  \code{\link[=coef.egf]{coef}(object, full = TRUE)}
  and preserved as an attribute of the result.}
\item{\dots}{
  unused optional arguments.}
}
\value{
A data frame with one row per coefficient and variables:
  
\item{cov}{
  label for a covariance matrix.
  This is the interaction of \code{term} and \code{group},
  but using levels with format \code{"Sigma[\%d]"}.}
\item{vec}{
  label for a random vector.
  This is the interaction of \code{term}, \code{group}, and \code{level},
  but using levels with format \code{"u[\%d]"}.}
\item{bottom}{
  label for a bottom level mixed effects model parameter,
  in this case for a random effect coefficient;
  this is a string with format \code{"b[\%d]"}.}
\item{top}{
  name of the top level nonlinear model parameter whose
  fitted value is a function of \code{bottom},
  from \code{\link{egf_top}(object)}.}
\item{term, group}{
  term from the random effects component of the mixed effects
  model formula for parameter \code{top}.
  \code{term} and \code{group} give the left and right hand sides
  of the term, which is a call to binary operator \code{|}.}
\item{level}{
  level of the factor or interaction indicated by \code{group}.}
\item{colname}{
  column name in the random effects design matrix
  \code{\link[=model.matrix.egf]{model.matrix}(object, "random")}.}
\item{value}{
  random effect conditional mode (unit variance scale),
  from segment \code{b} of
  \code{\link[=coef.egf]{coef}(object, full = TRUE)}.}

If \code{makeSigma = TRUE}, then the result has attribute \code{Sigma},
a list of covariance matrices corresponding to the levels of variable
\code{cov}.
}
\seealso{
The generic function \code{\link[nlme]{ranef}}.
}
