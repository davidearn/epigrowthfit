\name{simulate.egf}
\alias{simulate.egf}
\title{Simulation and parametric bootstrapping}
\usage{
\method{simulate}{egf}(
  object,
  nsim = 1L,
  seed = NULL,
  boot = FALSE,
  control = list(),
  parallel = egf_parallel(),
  trace = FALSE,
  ...
)
}
\arguments{
\item{object}{An \code{"\link{egf}"} object.}

\item{nsim}{A positive integer indicating a number of replications.}

\item{seed}{An integer used to set the \link{RNG} state before simulation.
The default (\code{NULL}) is to use the state at the time of
the function call. The RNG state (either a list of arguments to
\code{\link{set.seed}} or a value of \code{\link{.Random.seed}})
is preserved as an attribute of the result.}

\item{boot}{A \link{logical} flag. If \code{TRUE}, then a bootstrapping
step is performed.}

\item{control}{A list of control parameters passed to \code{\link{nlminb}}.}

\item{parallel}{An \code{"\link{egf_parallel}"} object defining options for \R level
parallelization.}

\item{trace}{A logical flag.
If \code{TRUE}, then basic tracing messages indicating progress
are printed.
Depending on \code{object$control$trace}, these may be mixed with
optimizer output.}

\item{...}{Unused optional arguments.}
}
\value{
A list inheriting from class \code{"egf_simulate"}, with elements:
\item{simulations}{
  A data frame containing simulated incidence data. It has variables
  \code{ts}, \code{window}, and \code{time}, and \code{nsim} further
  variables with names of the form \code{x[0-9]+}. It corresponds
  rowwise to \code{\link[=model.frame.egf]{model.frame}(object)}.
}
\item{boot}{
  If \code{boot = TRUE}, then a numeric matrix with \code{nsim} columns,
  each a sample from the conditional distribution of the parameter vector
  \code{unlist(\link[=coef.egf]{coef}(object, full = FALSE))}.
  Otherwise, \code{NULL}.
}
Attribute \code{RNGstate} preserves the RNG state prior to simulation,
making the result reproducible.
}
\description{
Simulates incidence data conditional on a fitted nonlinear
mixed effects model of epidemic growth.
(Only observations within fitting windows are simulated.)
Optionally re-estimates the model given the simulated data,
thus generating samples from the conditional distribution
of the bottom level parameter vector.
}
\details{
Bootstrap optimizations are typically expensive for nontrivial models.
They are parallelized at the C++ level when there is OpenMP support and
\code{object$control$omp_num_threads} is set to an integer greater than 1.
If there is no OpenMP support, then bootstrap optimizations can still be
parallelized at the \R level with appropriate setting of \code{parallel}.

Arguments \code{control}, \code{parallel}, and \code{trace} are unused
when \code{boot = FALSE}.
}
\examples{
object <- egf_cache("egf-1.rds")
zz <- egf_cache("simulate-egf-1.rds", {
    simulate(object, nsim = 6L, seed = 181952L, boot = TRUE)
})
str(zz)
matplot(t(zz$boot[object$nonrandom, ]), type = "o", las = 1,
        xlab = "simulation",
        ylab = "parameter value")

}
