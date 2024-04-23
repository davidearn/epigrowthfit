\name{getCall.egf}
\alias{getCall.egf}
\alias{getCall.egf_no_fit}
\title{Extract Model Calls}
\description{
Extracts from a model object the call to \code{\link{egf}} that produced it.
This method exists mainly to enable compatibility with the default method of
generic function \code{\link{update}}.
}
\usage{
\method{getCall}{egf}(x, \dots)
}
\arguments{
\item{x}{
  an \code{\link{egf}} object.}
\item{\dots}{
  unused optional arguments.}
}
\value{
A call to \code{\link{egf}}.
}
\seealso{
The generic function \code{\link{getCall}}.
}
