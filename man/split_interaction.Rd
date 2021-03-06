% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_formula.R
\name{split_interaction}
\alias{split_interaction}
\title{Split an interaction}
\usage{
split_interaction(x)
}
\arguments{
\item{x}{A \link{call}, \link{name}, or \link{atomic} scalar.}
}
\value{
A \link{list} of \link{call}s, \link{name}s, and \link{atomic} scalars.
}
\description{
Recursively constructs a \link{list} of arguments to binary operator
\code{`:`} from a nested \link{call} to \code{`:`}, excluding arguments
that are themselves calls to \code{`:`}.
}
\examples{
## x <- quote(a:b:I(f:g):log(h))
## split_interaction(x)

}
\keyword{internal}
