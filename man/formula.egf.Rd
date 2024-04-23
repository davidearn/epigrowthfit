\name{formula.egf}
\alias{formula.egf}
\alias{formula.egf_no_fit}
\title{Extract Model Formulae}
\description{
Extracts from a model object the mixed effects model formula
corresponding to a top level nonlinear model parameter.
}
\usage{
\method{formula}{egf}(x, top = egf_top(x), split = FALSE, \dots)
}
\arguments{
\item{x}{
  an \code{\link{egf}} object.}
\item{top}{
  a character string specifying a top level nonlinear model parameter.}
\item{split}{
  a logical flag.  If \code{TRUE}, then random effect terms are deleted
  from the formula and preserved as an attribute.}
\item{\dots}{
  unused optional arguments.}
}
\value{
By default, the mixed effects model formula corresponding to
\code{top}.  If \code{split = TRUE}, then the same formula
without random effect terms.  The deleted terms are stored in
an expression vector and preserved as attribute \code{random}
of the result.
}
\seealso{
The generic function \code{\link{formula}}.
}
