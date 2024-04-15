\name{simulate.egf}
\alias{simulate.egf}
\title{Simulation and Parametric Bootstrapping}
\description{
Simulates incidence data conditional on a fitted nonlinear mixed effects
model of epidemic growth.  Optionally re-estimates the model given the
simulated data, thus generating samples from the conditional distribution
of the bottom level parameter vector.
}
\usage{
\method{simulate}{egf}(object, nsim = 1, seed = NULL,
         bootstrap = FALSE,
         control = list(), parallel = egf_parallel(),
         trace = FALSE, \dots)
}
\arguments{
\item{object}{an \code{\link{egf}} object.}
\item{nsim}{a positive integer indicating a number of replications.}
\item{seed}{an integer used to set the \link{RNG} state before
  simulation or, otherwise, \code{NULL}; see \code{\link{simulate}}.}
\item{bootstrap}{a logical.  If \code{TRUE}, then a bootstrapping
  step is performed.}
\item{control}{passed to \code{\link{nlminb}}.}
\item{parallel}{an \code{\link{egf_parallel}} object defining options
  for \R{} level parallelization.}
\item{trace}{a logical.  If \code{TRUE}, then basic tracing messages
  indicating progress are printed.  These may be mixed with optimizer
  output depending on \code{object[["control"]][["trace"]]}.}
\item{\dots}{additional arguments passed from or to other methods.}
}
\value{
A list inheriting from class \code{simulate.egf}, with elements:

\item{simulation}{a data frame containing simulated incidence data.
  It has variables \code{ts}, \code{window}, \code{time}, and \code{X},
  where \code{X} is a numeric matrix with \code{nsim} columns.
  It corresponds rowwise to
  \code{\link[=model.frame.egf]{model.frame}(object)}.}
\item{bootstrap}{a numeric matrix with \code{nsim} columns, each a
  sample from the conditional distribution of the parameter vector
  represented by \code{\link[=coef.egf]{coef}(object)}.
  \code{NULL} if the method call did not set \code{bootstrap = TRUE}.}

Attribute \code{RNGstate} preserves the RNG state prior to simulation,
making the result reproducible.
}
\details{
Bootstrap optimizations are typically expensive for nontrivial models.
They are parallelized at the \proglang{C++} level when there is OpenMP
support and \code{object[["control"]][["omp_num_threads"]]} is set
to an integer greater than 1.  When there is no OpenMP support, they
can still be parallelized at the \R{} level with appropriate setting
of argument \code{parallel}.

Arguments \code{control}, \code{parallel}, and \code{trace} are unused
when \code{bootstrap = FALSE}.
}
\seealso{
\code{\link{simulate}}.
}
\examples{
object <- egf_cache("egf-2.rds")
zz <- egf_cache("simulate-egf-1.rds",
    simulate(object, nsim = 6L, seed = 181952L, bootstrap = TRUE))
str(zz)
matplot(t(zz[["bootstrap"]][!object[["random"]], ]),
        type = "o", las = 1, xlab = "simulation", ylab = "value")
}
