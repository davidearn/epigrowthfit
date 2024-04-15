\name{profile.egf}
\alias{profile.egf}
\alias{confint.profile.egf}
\alias{plot.profile.egf}
\title{Compute Univariate Likelihood Profiles}
\description{
Computes univariate likelihood profiles of fixed effect coefficients,
random effect covariance parameters, and linear combinations thereof.
}
\usage{
\method{profile}{egf}(fitted, A, which, level = 0.95, grid = 12L,
        parallel = egf_parallel(), trace = FALSE,
        top = egf_top(fitted), subset = NULL, select = NULL, \dots)

\method{confint}{profile.egf}(object, parm = seq_along(object), level = attr(object, "level"),
        maybeClass = FALSE, \dots)

\method{plot}{profile.egf}(x, parm = seq_along(x), level = attr(x, "level"),
     type = c("z", "abs(z)", "z^2"), \dots)
}
\arguments{
\item{fitted}{an \code{\link{egf}} object.}
\item{A}{a numeric matrix with \code{length(coef(fitted))} columns.
  Each row specifies a linear combination of the elements of
  \code{coef(fitted)} to be profiled.}
\item{which}{a valid index vector for \code{coef(fitted)}.
  Only indexed elements are profiled.  Setting \code{which = i}
  is equivalent to setting \code{A = diag(p)[i, ]} where
  \code{p = length(coef(fitted))}.}
\item{level}{a number in the interval \eqn{(0,1)} indicating a
  confidence level.  Profiles are computed up to a change in deviance
  equal to \code{qchisq(level, df = 1)}.}
\item{grid}{a positive integer.  Step sizes chosen adaptively by
  \code{\link[TMB]{tmbprofile}} will generate approximately this
  many points on each side of a profile's minimum point.}
\item{parallel}{an \code{\link{egf_parallel}} object defining options
  for \R{} level parallelization.}
\item{trace}{a logical.  If \code{TRUE}, then basic tracing messages
  indicating progress are printed.  These may be mixed with optimizer
  output depending on \code{fitted[["control"]][["trace"]]}.}
\item{top}{a subset of \code{\link{egf_top}(fitted)} naming top
  level nonlinear model parameters whose population fitted values
  (see \code{\link{fitted.egf}}) should be profiled.}
\item{subset, select}{index vectors for the rows and columns of
  \code{model.frame(fitted, "combined")} or language objects
  evaluating to such vectors.
  \code{subset} indicates fitting windows for which profiles
  should be computed; the default indicates all.
  \code{select} indicates variables that should be appended to
  the result; the default indicates none.
  Evaluation of language objects follows the implementation of
  \code{\link{subset.data.frame}}.}
\item{\dots}{additional arguments passed from or to other methods.}
\item{object, x}{a \code{profile.egf} object.}
\item{parm}{a valid index vector for \code{object} or \code{x}
  indicating a subset of the profiles.}
\item{maybeClass}{a logical.  If \code{TRUE} and
  if \code{object} was created by a call to \code{\link{profile}}
  not setting \code{A} or \code{which}, then the value of the
  \code{\link{confint}} call is a \code{\link{confint.egf}} object,
  not a matrix.}
\item{type}{a character string indicating which of \eqn{z},
  \eqn{|z|}{abs(z)}, and \eqn{z^{2}}{z^2} is plotted.}
}
\value{
A list of length \code{nrow(A)} inheriting from classes
\code{profile.egf} and \code{profile}.  Each element is a data frame
specifying a profile, with two variables:

\item{z}{a numeric vector containing profile \eqn{z}-statistics.
  The profile \eqn{z}-statistic is the appropriately signed square
  root of the change in deviance under the restricted model.}
\item{par.vals}{a numeric matrix with one column containing values
  of the linear combination specified by \code{A[i, ]}.}

The coefficient matrix \code{A} and confidence level \code{level}
are preserved as attributes.
}
\details{
Computation of likelihood profiles is typically expensive as it requires
estimation of many restricted models.
It is parallelized at the \proglang{C++} level when there is OpenMP
support and \code{fitted[["control"]][["omp_num_threads"]]} is set
to an integer greater than 1.  When there is no OpenMP support, it
can still be parallelized at the \R{} level with appropriate setting
of argument \code{parallel}.

Setting more than one of the arguments \code{A}, \code{which}, and
\code{top} in the method call is an error.
}
\seealso{
\code{\link{profile}}, \code{\link{plot.profile}}.
}
\examples{
fitted <- egf_cache("egf-1.rds")
zz <- egf_cache("profile-egf-1.rds",
    profile(fitted, top = "log(r)",
            subset = quote(country == "A" & wave == 1)))
str(zz)

confint(object)

pty <- c("z", "abs(z)", "z^2")
bty <- c("o", "u", "u")
for (i in 1:3)
    plot(object, type = pty[i], bty = bty[i], las = 1)
}
