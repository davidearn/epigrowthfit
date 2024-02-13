\name{plot.confint.egf}
\alias{plot.confint.egf}
\title{Plot confidence intervals}
\usage{
\method{plot}{confint.egf}(
  x,
  per_plot = 12L,
  subset = NULL,
  order = NULL,
  label = NULL,
  main = NULL,
  ...
)
}
\arguments{
\item{x}{An \code{"\link{confint.egf}"} object.}

\item{per_plot}{A positive integer. One plot will display at most this many
confidence intervals.}

\item{subset}{An index vector for \code{seq_len(nrow(x))}
or a language object evaluating in \code{x} to such a vector.
\code{subset} indicates which confidence intervals are plotted;
the default (\code{NULL}) indicates all of them.}

\item{order}{A permutation of \code{seq_len(nrow(x))}
or a language object evaluating in \code{x} to such a vector.
\code{order} indicates the order in which confidence intervals are plotted;
the default (\code{NULL}) indicates the original order.}

\item{label}{A character or expression vector of length \code{nrow(x)}
(or a language object evaluating in \code{x} to such a vector).
\code{label} indicates \eqn{y}-axis labels for confidence intervals;
the default (\code{NULL}) is to obtain labels from \code{x[["window"]]}.}

\item{main}{A character or expression vector of length 1
indicating a plot title, to be recycled for all plots.}

\item{...}{Unused optional arguments.}
}
\value{
\code{NULL} (invisibly).
}
\description{
A method for graphically comparing confidence intervals
on fitted values of top level nonlinear model parameters.
It creates a one-dimensional plot with
confidence intervals drawn as stacked horizontal line segments.
}
\examples{
x <- egf_cache("confint-egf-1.rds")

op <- par(mar = c(4.5, 4, 2, 1), oma = c(0, 0, 0, 0))
plot(x)
par(op)
}
