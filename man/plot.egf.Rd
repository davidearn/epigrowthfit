\name{plot.egf}
\alias{plot.egf}
\title{Plot nonlinear mixed effects models of epidemic growth}
\usage{
\method{plot}{egf}(
  x,
  type = c("interval", "cumulative", "rt"),
  time_as = c("Date", "numeric"),
  dt = 1,
  log = TRUE,
  zero = NA,
  show_predict = TRUE,
  show_tdoubling = FALSE,
  level = 0.95,
  control = egf_control_plot(),
  cache = NULL,
  plot = TRUE,
  subset = NULL,
  order = NULL,
  xlim = NULL,
  ylim = NULL,
  main = NULL,
  sub = NULL,
  xlab = NULL,
  ylab = NULL,
  ylab_outer = NULL,
  ...
)
}
\arguments{
\item{x}{An \code{"\link{egf}"} object.}

\item{type}{A character string indicating a type of plot. The options are:
interval incidence (\code{"interval"}),
cumulative incidence (\code{"cumulative"}), and
per capita growth rate (\code{"rt"}).}

\item{time_as}{A character string indicating how numeric times are displayed
on the bottom axis. The options are:
as is (\code{"numeric"})
and with a calendar (\code{"Date"}).
In the latter case, horizontal user coordinates on measure time in days
since \code{1970-01-01 00:00:00}.}

\item{dt}{A positive number indicating an observation interval.
Predicted curves are evaluated on grids with this spacing.
When \code{type = "interval"}, counts observed over shorter
or longer intervals \code{dt0} are scaled by a factor of
\code{dt / dt0} so that their scale matches that of the curves.
Scaled counts can be highlighted via \code{control}.
If \code{x} specifies a model with day of week effects
(\code{x$model$day_of_week > 0}), then setting \code{dt}
has no effect as it is set to 1 internally.}

\item{log}{A logical flag. If \code{TRUE}, then the dependent variable
is plotted on a logarithmic scale.}

\item{zero}{A positive number indicating a line on which to plot zeros
when \code{log = TRUE}.
\code{NA} is to place zeros on the bottom axis.
\code{NULL} is to suppress zeros.
\code{zero} is used for \code{type = "interval"} or \code{"cumulative"} only.}

\item{show_predict}{An integer flag:
2 is to draw predicted curves with confidence bands,
1 is draw predicted curves only,
0 is to draw neither.
Logical values are coerced to integer.}

\item{show_tdoubling}{An integer flag:
2 is to print initial doubling time estimates in the top margin
with confidence intervals,
1 is to print estimates only,
0 is to print neither.
Logical values are coerced to integer.
Supported only if \code{x$model$curve} is
\code{"exponential"}, \code{"logistic"}, or \code{"richards"}.}

\item{level}{A number in the interval (0,1). This is the confidence level used
when \code{show_predict = 2} or \code{show_tdoubling = 2}.}

\item{control}{An \code{"\link{egf_control_plot}"} object controlling the appearance
of most plot elements.}

\item{cache}{A \code{"plot.egf"} object returned by a previous evaluation
of \code{plot(x)}. Fitted and predicted values and standard errors
stored in \code{cache} are reused if possible to avoid recomputation.}

\item{plot}{A \link{logical} flag. If \code{FALSE}, then nothing is plotted.
Useful when only the returned \code{"plot.egf"} object is desired.}

\item{subset}{An index vector for the rows of
\code{mf = \link[=model.frame.egf]{model.frame}(object, "combined")}
or a language object evaluating to such a vector.
Only time series corresponding to indexed rows are plotted and
only fitting windows corresponding to indexed rows are highlighted;
the default (\code{NULL}) is to plot all time series and highlight
all fitting windows.
Evaluation of language objects follows the implementation of
\code{\link{subset}(mf, subset, select)}.}

\item{order}{A permutation of \code{seq_len(nrow(mf))}
or a language object evaluating in \code{mf} to such a vector.
\code{order} indicates the order in which time series are plotted;
the default (\code{NULL}) indicates the original order.}

\item{xlim, ylim}{Numeric vectors of length 2 specifying axis limits,
which are recycled for all plots.
If \code{time_as = "Date"}, then \code{xlim} can instead
be a \link{Date} vector or any vector coercible to Date
via \code{\link{as.Date}(xlim)}.}

\item{main, sub, xlab, ylab, ylab_outer}{Character or expression
vectors or (\code{main}, \code{sub}) language objects evaluating
in \code{mf} to such vectors.
These are used to generate plot (\code{main}, \code{sub}) and
axis (\code{xlab}, \code{ylab}, \code{ylab_outer}) labels.
\code{ylab_outer} is used for \code{type = "rt"} only.}

\item{...}{Unused optional arguments.}
}
\value{
A data.frame inheriting from class \code{"plot.egf"}.
If argument \code{cache} was supplied in the function call,
then this data frame is the result of augmenting \code{cache}
with new computations.
}
\description{
Methods for plotting \code{"\link{egf}"} objects.
}
\details{
Computation of fitted and predicted values and standard errors
is performed before any plots are created. To avoid waste of
computation time, cached computations are returned \emph{even if}
an error is thrown during plotting. To ensure that the cache is
preserved, assign the result of the function call to a name:
\code{cache <- plot(x, \dots)}.

Caching functionality must be used with care, as mismatch between
\code{x} and \code{cache} will not be detected. Constructions such
as \code{plot(y, cache = plot(x, \dots), \dots)}, where \code{x}
and \code{y} are different objects, should not be expected to produce
correct results.
}
\examples{
x <- egf_cache("egf-1.rds")
l <- list(legend = list(cex = 0.8),
          estimate = list(cex = 0.8, font = 2),
          ci = list(cex = 0.8))
control <- egf_control_plot(tdoubling = l)

op <- par(mar = c(3.5, 5, 5, 1))
plot(x,
     type = "interval",
     show_predict = 2L,
     show_tdoubling = 2L,
     control = control)
plot(x,
     type = "cumulative",
     main = "Fitted exponential model",
     sub = quote(paste("Country", country)))
par(op)

op <- par(mar = c(3.5, 9.5, 5, 1))
plot(x, type = "rt", subset = quote(country \%in\% LETTERS[4:6]))
par(op)

}
