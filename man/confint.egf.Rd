\name{confint.egf}
\alias{confint.egf}
\alias{plot.confint.egf}
\title{Confidence Intervals}
\description{
Computes confidence intervals on fixed effect coefficients,
random effect covariance parameters, and linear combinations thereof,
including population fitted values.
Intervals on individual fitted values accounting for random effects
are supported, but only by \code{method = "wald"}.
}
\usage{
\method{confint}{egf}(object, parm, level = 0.95,
        A = seq_along(par), method = c("wald", "profile", "uniroot"), scale = 7,
        parallel = egf_parallel(), trace = FALSE,
        top = egf_top(object), subset = NULL, select = NULL,
        class = FALSE, link = TRUE, random = FALSE, \dots)

\method{plot}{confint.egf}(x, by = 12L,
     subset = NULL, order = NULL, label = NULL, main = NULL, \dots)
}
\arguments{
\item{object}{
  an \code{\link{egf}} object.}
\item{parm}{
  unused argument, for consistency with the generic function.}
\item{level}{
  a number in the interval \eqn{(0,1)} indicating a confidence level.}
\item{A}{
  a numeric matrix with \code{1+p} columns,
  where \code{p = length(coef(object))},
  in which case each row specifies a linear combination of the elements
  of \code{c(1, coef(fitted))} for which intervals are computed;
  or a valid index vector for \code{coef(fitted)},
  in which case intervals are computed for the indexed elements;
  or \code{NULL},
  in which case intervals on fitted values are computed.}
\item{method}{
  a character string indicating \emph{how} intervals are computed.}
\item{scale}{
  a positive number, for \code{method = "uniroot"}.
  \code{\link[TMB]{tmbroot}} will search for roots between
  \code{value-scale*se} and \code{value+scale*se}, where
  \code{value} and \code{se} are the estimate and standard error.}
\item{parallel}{
  an \code{\link{egf_parallel}} object defining options for \R{} level
  parallelization.}
\item{trace}{
  a logical.  If \code{TRUE}, then basic tracing messages indicating
  progress are printed.  These may be mixed with optimizer output
  depending on \code{object[["control"]][["trace"]]}.}
\item{top}{
  a subset of \code{\link{egf_top}(object)} naming top level
  nonlinear model parameters for which intervals on fitted values
  should be computed.}
\item{subset, select}{
  index vectors for the rows and columns of
  \code{model.frame(object, "combined")} or language objects
  evaluating to such vectors.
  \code{subset} indicates fitting windows for which intervals
  should be computed; the default indicates all.
  \code{select} indicates variables that should be appended to
  the result; the default indicates none.
  Evaluation of language objects follows the implementation of
  \code{\link{subset.data.frame}}.\cr
  \cr
  For the \code{\link{plot}} method:\cr
  an index vector for \code{seq_len(nrow(x))} or a language object
  evaluating in \code{x} to such a vector.
  \code{subset} indicates which intervals are plotted;
  the default indicates all.}
\item{class}{
  a logical.  If \code{TRUE} and if \code{A = NULL},
  then the value of the \code{\link{confint}} call is
  a \code{\link{confint.egf}} object, not a matrix.}
\item{link}{
  a logical. If \code{FALSE} and if \code{A = NULL}
  and \code{class = TRUE},
  then fitted values and confidence limits are returned on the
  \dQuote{natural} (inverse link) scale.}
\item{random}{
  a logical, affecting only \code{method = "wald"}.  If \code{TRUE},
  then intervals are computed for individual fitted values,
  which count random effects, rather than population fitted values,
  which do not.}
\item{\dots}{
  additional arguments passed from or to other methods.}
\item{x}{
  a \code{\link{confint.egf}} object.}
\item{by}{
  a positive integer indicating the number of intervals displayed
  in one plot.}
\item{order}{
  a permutation of \code{seq_len(nrow(x))} or a language object
  evaluating in \code{x} to such a vector.
  \code{order} indicates the order in which intervals are plotted;
  the default indicates the original order.}
\item{label}{
  a character or expression vector of length \code{nrow(x)}
  or a language object evaluating in \code{x} to such a vector.
  \code{label} indicates \eqn{y}-axis labels for intervals;
  the default is to use \code{as.character(x[["window"]])}.}
\item{main}{
  a character or expression vector of length 1 indicating a plot title,
  to be recycled for all plots.}
}
\value{
A numeric array in 2 or 3 dimensions containing the lower and upper
confidence limits in the last dimension.

When confidence intervals on fitted values are desired, the user will
set \code{A = NULL} and in that case have the option of passing
\code{class = TRUE} to obtain an augmented result.  Thus, alternatively:

A data frame inheriting from class \code{confint.egf}, with variables:

\item{top}{
  top level nonlinear model parameter, from
  \code{\link{egf_top}(object)}.}
\item{ts}{
  time series, from
  \code{levels(\link[=model.frame.egf]{model.frame}(object)[["ts"]])}.}
\item{window}{
  fitting window, from
  \code{levels(\link[=model.frame.egf]{model.frame}(object)[["window"]])}.}
\item{value}{
  fitted value.}
\item{ci}{
  a numeric matrix with two columns giving the lower and upper
  confidence limits.}
\item{\dots}{
  further variables from
  \link[=model.frame.egf]{model.frame}(object, "combined")
  specified by argument \code{select}.}

The confidence level \code{level} is preserved as an attribute.
}
\details{
Three methods for computing confidence intervals are available:
\describe{
\item{\code{"wald"}}{
  confidence limits are calculated as
  \preformatted{value + c(-1, 1) * sqrt(q) * se}
  where \code{q = qchisq(level, df = 1)}.}
\item{\code{"profile", "uniroot"}}{
  confidence limits are calculated as approximate solutions of the
  equation
  \preformatted{2 * (f(x) - f(value)) = q}
  where \code{q = qchisq(level, df = 1)}
  and \code{f} is the negative log marginal likelihood function
  expressed as a function of the parameter \code{x} in question.
  Solutions are approximated by interpolating a likelihood profile
  (\code{"profile"}) or by rootfinding (\code{"uniroot"}).}
}
\code{"wald"} assumes asymptotic normality of the maximum likelihood
estimator.  \code{"profile"} and \code{"uniroot"} avoid this contraint
but are typically expensive, requiring estimation of many restricted
models.
They are parallelized at the \proglang{C++} level when there is OpenMP
support and \code{object[["control"]][["omp_num_threads"]]} is set
to an integer greater than 1.  When there is no OpenMP support, they
can still be parallelized at the \R{} level with appropriate setting
of argument \code{parallel}.
}
\seealso{
The generic function \code{\link{confint}}.
}
\examples{
object <- egf_cache("egf-1.rds")
zz1 <- egf_cache("confint-egf-1.rds",
    confint(object, A = NULL, method = "wald", class = TRUE,
            random = TRUE))
str(zz1)

op <- par(mar = c(4.5, 4, 2, 1), oma = c(0, 0, 0, 0))
plot(zz1)
par(op)

\donttest{
zz2 <- egf_cache("confint-egf-2.rds",
    confint(object, A = NULL, method = "profile", class = TRUE,
            top = "log(r)", subset = quote(country == "A" & wave == 1)))
zz3 <- egf_cache("confint-egf-3.rds",
    confint(object, A = NULL, method = "uniroot", class = TRUE,
            top = "log(r)", subset = quote(country == "A" & wave == 1)))
}
}
