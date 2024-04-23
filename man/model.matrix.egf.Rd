\name{model.matrix.egf}
\alias{model.matrix.egf}
\alias{model.matrix.egf_no_fit}
\title{Extract Design Matrices}
\description{
Extracts from a model object fixed and random effects design matrices.
}
\usage{
\method{model.matrix}{egf}(object, which = c("fixed", "random"),
             top = NULL, random = NULL, \dots)
}
\arguments{
\item{object}{
  an \code{\link{egf}} object.}
\item{which}{
  a character string controlling what is returned:
  \describe{
  \item{\code{"fixed"}}{
    the fixed effects design matrix \code{X} corresponding to \code{top}.}
  \item{\code{"random"}}{
    the random effects design matrix \code{Z} corresponding to \code{top}
    and \code{random}.}
  }}
\item{top}{
  a character string specifying a top level nonlinear model parameter.
  \code{NULL} indicates \dQuote{all of them}; see \sQuote{Details}.}
\item{random}{
  a random effect term, which is a call to binary operator \code{|}.
  \code{NULL} indicates \dQuote{all of them}; see \sQuote{Details}.}
\item{\dots}{
  unused optional arguments.}
}
\value{
A (sparse) \link[Matrix:dgCMatrix-class]{dgCMatrix} or a traditional
(dense) matrix, with attributes \code{assign} and \code{contrasts}
except in special cases; see \sQuote{Details}.
}
\details{
\code{model.matrix(which = "fixed", top = NULL)}
returns the result of combining (in the sense of \code{cbind})
all fixed effects design matrices.

\code{model.matrix(which = "random", top = "<name>", random = NULL)}
returns the result of combining all random effects design matrices
associated with parameter \code{top}.

\code{model.matrix(which = "random", top = NULL, random = NULL)}
returns the result of combining all random effects design matrices
\emph{and} permuting the columns to obtain a convenient ordering of
random effect coefficients.
(Coefficients are sorted by relation to a common random vector.
Random vectors are sorted by relation to a common covariance matrix.)

None of these \dQuote{combined} design matrices possesses attributes
\code{assign} and \code{contrasts}.
}
\seealso{
The generic function \code{\link{model.matrix}}.
}
