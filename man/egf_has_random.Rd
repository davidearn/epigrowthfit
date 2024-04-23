\name{egf_has_random}
\alias{egf_has_random}
\title{Test for Random Effects}
\description{
Tests whether an object specifies a model with random effects.
}
\usage{
egf_has_random(object, check = TRUE)
}
\arguments{
\item{object}{
  an \code{\link{egf}} object.}
\item{check}{
  a logical.  If \code{TRUE}, then an error is thrown if \code{object}
  does not actually inherit from class \code{\link{egf}}.}
}
\value{
\code{TRUE} or \code{FALSE}.
}
