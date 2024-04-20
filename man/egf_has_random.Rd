\name{egf_has_random}
\alias{egf_has_random}
\title{Check for random effects}
\usage{
egf_has_random(object, check = TRUE)
}
\arguments{
\item{object}{An \code{"\link{egf}"} or \code{"egf_no_fit"} object.}
\item{check}{a logical.  If \code{TRUE}, then test inheritance.}
}
\value{
\code{TRUE} or \code{FALSE}.
}
\description{
Determines whether an object specifies a random effects model.
}
