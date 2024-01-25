\name{plot.egf_confint}
\alias{plot.egf_confint}
\title{Plot confidence intervals}
\usage{
\method{plot}{egf_confint}(
  x,
  per_plot = 12L,
  subset = NULL,
  order = NULL,
  main = NULL,
  label = NULL,
  ...
)
}
\arguments{
\item{x}{An \code{"\link[=confint.egf]{egf_confint}"} object.}

\item{per_plot}{A positive integer. One plot will display at most this many
confidence intervals.}

\item{subset}{An expression to be evaluated in \code{x}.
It must evaluate to a valid index vector for the rows of \code{x}.
Only indexed confidence intervals are plotted.
The default (\code{NULL}) is to plot all confidence intervals.}

\item{order}{An expression to be evaluated in \code{x},
typically a call to \code{\link{order}},
determining the order in which confidence intervals
are plotted.
It must evaluate to a permutation of \code{seq_len(nrow(x))}.
The default (\code{NULL}) is equivalent to \code{seq_len(nrow(x))}.}

\item{main}{An expression or character string indicating a plot title,
to be recycled for all plots.}

\item{label}{An expression to be evaluated in \code{x}, typically a factor,
an interaction of factors, or a call to \code{\link{sprintf}}.
It is used to create appropriate \eqn{y}-axis labels for confidence
intervals. The default (\code{NULL}) is to take labels from \code{x$window}.}

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
\details{
See topic \code{\link{egf_eval}} for details on nonstandard evaluation
of \code{subset}, \code{order}, and \code{label}.
}
\examples{
x <- egf_cache("confint-egf-1.rds")

op <- par(mar = c(4.5, 4, 2, 1), oma = c(0, 0, 0, 0))
plot(x)
par(op)
}
