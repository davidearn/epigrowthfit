\name{egf_control}
\alias{egf_control}
\title{Define Control Parameters}
\description{
Set parameters controlling the behaviour of \code{\link{egf}}.
}
\usage{
egf_control(outer_optimizer = egf_optimizer(nlminb),
            inner_optimizer = egf_optimizer(newton),
            trace = FALSE, profile = FALSE, sparse_X = FALSE,
            omp_num_threads = getOption("egf.cores", 1L))
}
\arguments{
\item{outer_optimizer, inner_optimizer}{
  \code{\link{egf_optimizer}} objects specifying \dQuote{outer} and
  \dQuote{inner} optimization methods.}
\item{trace}{
  an integer determining the amount of tracing performed;
  see \sQuote{Details}.}
\item{profile}{
  a logical.  If \code{TRUE}, then fixed effect coefficients are
  profiled out of the likelihood, which may stabilize optimization
  for models with many fixed effects.
  This \dQuote{feature} should be considered experimental,
  and in fact it may \emph{de}stabilize optimization, as it relies
  on assumptions about the optimization problem that are not
  necessarily satisfied by the nonlinear mixed effects models fit
  by \code{\link{egf}}.}
\item{sparse_X}{
  a logical.  If \code{TRUE}, then the fixed effects design matrix
  is represented as a (sparse) \link[Matrix:dgCMatrix-class]{dgCMatrix},
  rather than as a traditional (dense) matrix.}
\item{omp_num_threads}{
  an integer indicating a number of OpenMP threads to be used
  when evaluating the objective function, provided that
  \pkg{epigrowthfit} was compiled with OpenMP support.}
}
\value{
A list inheriting from class \code{egf_control} containing the
validated arguments.
}
\details{
\code{trace} affects the amount of information printed during
likelihood evaluations:

\describe{
\item{\code{0}}{
  likelihood evaluations are always silent.}
\item{\code{1}}{
  a message is printed whenever a negative log marginal likelihood
  term is \code{NaN} or exceeds \code{1e+09}.}
\item{\code{2}}{
  all negative log marginal likelihood terms are printed.}
}

\code{\link{egf}} passes \code{silent = trace == 0L} to
\code{\link[TMB]{MakeADFun}}.  A corollary is that nonzero values of
\code{trace} have a number of additional side effects:

\itemize{
\item error messages are printed during function and gradient evaluations;
\item the maximum absolute gradient element is printed with each gradient
evaluation; and
\item trace flags set by \code{\link[TMB]{config}} are turned on.
}
}
\section{Warning}{
Setting \code{trace > 0L} and \code{omp_num_threads > 0L} simultaneously
should be avoided, because tracing messages are printed using \R{} API
functions that are not thread-safe.
}
\examples{
control <- egf_control()
str(control)
}
