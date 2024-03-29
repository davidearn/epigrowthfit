\name{egf_optimizer}
\alias{egf_optimizer}
\alias{egf_inner_optimizer}
\title{Define an optimization method}
\usage{
egf_optimizer(f = nlminb, args = list(), control = list())

egf_inner_optimizer(f = newton, args = list(), control = list())
}
\arguments{
\item{f}{A function performing optimization. The outer optimization permits
\code{\link{optim}}, \code{\link{nlminb}}, and \code{\link{nlm}}
and any \code{optim}-like function. An \code{optim}-like function
is a function \code{f} such that:
(i) the first three arguments of \code{f} specify an initial parameter
vector, an objective function, and a gradient function, respectively;
(ii) \code{f} accepts \code{control} as a fourth (or later) argument;
and
(iii) \code{f} returns a list with elements \code{par}, \code{value},
\code{convergence}, and \code{message}.
The inner optimization permits \code{optim} and \code{\link[TMB]{newton}}
only.}

\item{args}{A list of arguments to \code{f} other than \code{control}.
If \code{f = \link{optim}} and \code{args} does not have \code{method}
as an element, then \code{method = "BFGS"} is appended.}

\item{control}{A list of control parameters to be assigned to argument \code{control}
of \code{f}.}
}
\value{
\code{egf_optimizer} returns a list inheriting from class
\code{"egf_optimizer"}, with elements:
\item{f}{
  An \code{\link{optim}}-like \link{function}, typically the result
  of wrapping the supplied optimizer to make it \code{optim}-like.
}
\item{args}{
  The supplied list of arguments
  (after possible deletion of elements with reserved names).
}
\item{control}{
  The supplied list of control parameters.
}

\code{egf_inner_optimizer} returns a list inheriting from class
\code{"egf_inner_optimizer"}, with elements:
\item{method}{
  A character string. This is \code{args$method} if \code{f = \link{optim}}
  and \code{"newton"} if \code{f = \link[TMB]{newton}}.
}
\item{control}{
  A list. This is \code{control} if \code{f = \link{optim}} and \code{args}
  (after possible deletion of elements with reserved names)
  if \code{f = \link[TMB]{newton}}.
  To align the default behaviour of \code{newton} with that of \code{optim},
  \code{trace = 0} is set if not specified in \code{args}.
}
}
\description{
These two functions link an optimizer with function arguments and
control parameters to define an optimization method for use by
\code{\link{egf}}. "Outer" and "inner" optimization methods are
defined separately by \code{egf_optimizer} and \code{egf_inner_optimizer}.
}
\examples{
optimizer <- egf_optimizer()
str(optimizer)

inner_optimizer <- egf_inner_optimizer()
str(inner_optimizer)

}
\seealso{
\code{\link[TMB]{MakeADFun}} for some details about outer and inner
optimizations

Other behaviour-defining functions: 
\code{\link{egf_control}()},
\code{\link{egf_model}()},
\code{\link{egf_parallel}()}

Other behaviour-defining functions: 
\code{\link{egf_control}()},
\code{\link{egf_model}()},
\code{\link{egf_parallel}()}
}
\concept{behaviour-defining functions}
