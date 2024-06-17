\name{fitted.egf}
\alias{fitted.egf}
\alias{fitted.egf_no_fit}
\alias{confint.fitted.egf}
\title{Fitted Values}
\description{
Retrieves fitted values of top level nonlinear model parameters.
The fitted value of a given parameter for a given fitting window
is obtained by adding
(i) the population fitted value computed as a linear combination
of fixed effect coefficients and
(ii) all applicable random effects, with random effects set equal
to their conditional modes.
}
\usage{
\method{fitted}{egf}(object,
       top = egf_top(object), subset = NULL, select = NULL,
       class = FALSE, se = FALSE, \dots)

\method{confint}{fitted.egf}(object, parm = seq_len(nrow(object)), level = 0.95,
        class = FALSE, \dots)
}
\arguments{
\item{object}{
  an \code{\link{egf}} or \code{fitted.egf} object.}
\item{top}{
  a subset of \code{\link{egf_top}(object)} naming top level
  nonlinear model parameters whose fitted values should be
  retrieved.}
\item{subset, select}{
  index vectors for the rows and columns of
  \code{model.frame(object, "combined")} or language objects
  evaluating to such vectors.
  \code{subset} indicates fitting windows for which fitted values
  should be retrieved; the default indicates all.
  \code{select} indicates variables that should be appended to
  the result; the default indicates none.
  Evaluation of language objects follows the implementation of
  \code{\link{subset.data.frame}}.}
\item{class}{
  a logical.  If \code{TRUE}, then the value of the method call
  is a \code{fitted.egf} or \code{\link{confint.egf}} object,
  not a matrix.}
\item{se}{
  a logical.  If \code{TRUE} and if \code{class = TRUE},
  then the result is augmented with approximate delta method
  standard errors.}
\item{\dots}{
  additional arguments passed from or to other methods.}
\item{parm}{
  a valid index vector for the rows of \code{object} indicating a subset
  of the fitted values.}
\item{level}{
  a number in the interval \eqn{(0,1)} indicating a confidence level.}
}
\value{
A numeric matrix containing fitted values.

Users can pass \code{class = TRUE} to obtain an augmented result.
Thus, alternatively:

A data frame inheriting from class \code{fitted.egf}, with variables:

\item{top}{
  top level nonlinear model parameter, from
  \code{\link{egf_top}(object)}.}
\item{ts}{
  time series, from
  \code{levels(\link[=model.frame.egf]{model.frame}(object)[["ts"]])}.}
\item{window}{
  fitting window, from
  \code{levels(\link[=model.frame.egf]{model.frame}(object)[["window"]])}.}
\item{value}{
  fitted value.}
\item{se}{
  approximate delta method standard error (only if requested).}
\item{\dots}{
  further variables from
  \link[=model.frame.egf]{model.frame}(object, "combined")
  specified by argument \code{select}.}
}
\seealso{
The generic function \code{\link{fitted}}.
}
\examples{
\donttest{
example("egf", package = "epigrowthfit")

zz <- fitted(m1, class = TRUE, se = TRUE)
str(zz)

confint(zz, class = TRUE)
}
}
