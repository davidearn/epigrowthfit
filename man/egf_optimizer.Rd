\name{egf_optimizer}
\alias{egf_optimizer}
\title{Define an Optimization Method}
\description{
Utilities for linking an optimizer with optional arguments and
control parameters to define an optimization method for use by
\code{\link{egf}}.
}
\usage{
egf_optimizer(f = nlminb, args = list(), control = list())
}
\arguments{
\item{f}{
  a function performing optimization.  Supported are
  \code{\link[TMB]{newton}}, \code{\link{optim}}, \code{\link{nlminb}},
  \code{\link{nlm}}, and any \code{optim}-like function.}
\item{args}{
  a list of optional arguments to \code{f} not including
  \code{control}.  If \code{f = \link{optim}} and \code{args} does
  not have \code{method} as an element, then \code{method = "BFGS"}
  is appended.}
\item{control}{
  a list of control parameters to be passed to \code{f}.}
}
\value{
A list inheriting from class \code{egf_optimizer} containing the
validated arguments, wherein \code{f} may be a new function wrapping
the supplied one to make it \code{\link{optim}}-like.
}
\details{
An \code{\link{optim}}-like function is a function \code{f} such that:

\itemize{  
\item the first three arguments of \code{f} specify an initial
parameter vector, an objective function, and a gradient function,
respectively;
\item \code{f} accepts \code{control} as a fourth (or later) argument;
and
\item \code{f} returns a named list with elements
\code{par}, \code{value}, \code{convergence}, and \code{message}.
}
}
\examples{
optimizer <- egf_optimizer(nlminb)
str(optimizer)
}
