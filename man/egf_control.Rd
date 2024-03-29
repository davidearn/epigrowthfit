\name{egf_control}
\alias{egf_control}
\title{Define control parameters}
\usage{
egf_control(
  optimizer = egf_optimizer(),
  inner_optimizer = egf_inner_optimizer(),
  trace = FALSE,
  profile = FALSE,
  sparse_X = FALSE,
  omp_num_threads = getOption("egf.cores", 1L)
)
}
\arguments{
\item{optimizer}{An \code{"\link{egf_optimizer}"} object, specifying an "outer"
optimization method.}

\item{inner_optimizer}{An \code{"\link{egf_inner_optimizer}"} object, specifying an "inner"
optimization method, or a list of such objects, in which case the
listed methods are tried in order until one succeeds. (If none succeed,
then a warning is issued.)}

\item{trace}{An integer flag determining the amount of tracing performed
(see Details). Logical values are coerced to integer.}

\item{profile}{A logical flag. If \code{TRUE}, then fixed effect parameters are profiled
out of the likelihood, which may stabilize optimization for models
with many fixed effects. This feature should be considered experimental,
and in fact may \emph{de}stabilize optimization, as it may rely on
assumptions about the optimization problem that are not necessarily
satisfied by the nonlinear mixed effects models fit by \code{\link{egf}}.}

\item{sparse_X}{A logical flag. If \code{TRUE}, then the fixed effects
\link[=model.matrix]{design} matrix is constructed in
\link[Matrix:sparseMatrix]{sparse} format.}

\item{omp_num_threads}{An integer specifying a number of OpenMP threads to be used
(if supported) when evaluating the objective function.}
}
\value{
A list inheriting from class \code{"egf_control"} containing the arguments
(after possible coercion).
}
\description{
Set parameters controlling the behaviour of \code{\link{egf}}.
}
\details{
\code{trace} affects the amount of information printed during
likelihood evaluations:
\describe{
\item{0}{
  Likelihood evaluations are always silent.
}
\item{1}{
  A message is printed whenever a negative log likelihood term
  is non-finite or exceeds \code{1e+09}.
}
\item{2}{
  All negative log likelihood terms are printed.
}
}

\code{\link{egf}} passes \code{silent = (trace == 0L)}
to \code{\link[TMB]{MakeADFun}}. As a result, nonzero values
of \code{trace} have a number of additional side effects:
\itemize{
\item error messages are printed during function and gradient evaluations;
\item the maximum absolute gradient element is printed with each gradient
evaluation; and
\item trace flags set by \code{\link[TMB]{config}} are turned on.
}
}
\section{Warning}{

Setting \code{trace > 0L} and \code{omp_num_threads > 0L}
simultaneously should be considered dangerous on builds of
\pkg{epigrowthfit} obtained from CRAN. These builds print
using R API that is not thread-safe. R API is avoided on
builds from source.
}

\examples{
control <- egf_control()
str(control)

}
\seealso{
Other behaviour-defining functions: 
\code{\link{egf_model}()},
\code{\link{egf_optimizer}()},
\code{\link{egf_parallel}()}
}
\concept{behaviour-defining functions}
