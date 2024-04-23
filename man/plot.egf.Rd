\name{plot.egf}
\alias{plot.egf}
\title{Plot Nonlinear Mixed Effects Models of Epidemic Growth}
\description{
A method for printing objects of class \code{\link{egf}}.
}
\usage{
\method{plot}{egf}(x, type = c("interval", "cumulative", "rt"),
     time_as = c("Date", "numeric"), delta = 1, log = TRUE, zero = NA,
     show_predict = TRUE, show_doubling = FALSE, level = 0.95,
     control = egf_control_plot(), cache = NULL, plot = TRUE,
     subset = NULL, order = NULL, xlim = NULL, ylim = NULL,
     main = NULL, sub = NULL, xlab = NULL, ylab = NULL, \dots)
}
\arguments{
\item{x}{
  an \code{\link{egf}} object.}
\item{type}{
  a character string indicating a type of plot.  The options are:
  interval incidence (\code{"interval"}),
  cumulative incidence (\code{"cumulative"}), and
  per capita growth rate (\code{"rt"}).}
\item{time_as}{
  a character string indicating how numeric times are displayed
  on the bottom axis.  The options are:
  as is (\code{"numeric"}) and
  with a calendar (\code{"Date"}).
  In the latter case, horizontal user coordinates on measure time
  in days since midnight on January 1, 1970.}
\item{delta}{
  a positive number indicating a step size on the time axis.
  Predicted curves are evaluated on a grid with this spacing.
  When \code{type = "interval"}, counts observed over shorter
  or longer intervals \code{delta0} are scaled by a factor
  of \code{delta/delta0} so that their scale matches that of
  predicted incidence.  Scaled counts can be highlighted via
  \code{control}.  If \code{x} specifies a model with day of
  week effects, then \code{delta = 1} is used unconditionally.}
\item{log}{
  a logical.  If \code{TRUE}, then the dependent variable is plotted
  on a logarithmic scale.}
\item{zero}{
  a positive number indicating a line on which to plot zeros when
  \code{log = TRUE} and \code{type = "interval"} or \code{"cumulative"}.
  \code{NA} is to place zeros on the bottom axis.
  \code{NULL} is to suppress zeros.}
\item{show_predict}{
  an integer flag:
  2 is to draw predicted curves with confidence bands,
  1 is draw predicted curves only,
  0 is to draw neither.}
\item{show_doubling}{an integer flag:
  2 is to print initial doubling time estimates in the top margin
  with confidence intervals,
  1 is to print estimates only,
  0 is to print neither.
  Nothing is printed for models without a well-defined initial
  exponential growth rate.}
\item{level}{
  a number in the interval \eqn{(0,1)} indicating confidence level,
  used when \code{show_predict = 2} or \code{show_doubling = 2}.}
\item{control}{
  an \code{\link{egf_control_plot}} object controlling the
  appearance of most plot elements.}
\item{cache}{
  a \code{plot.egf} object returned by a previous evaluation of
  \code{plot(x, \dots)}.  Fitted and predicted values and standard
  errors stored in \code{cache} are reused to avoid recomputation.}
\item{plot}{
  a logical.  If \code{FALSE}, then plotting does not occur.
  Useful when only the returned \code{plot.egf} object is desired.}
\item{subset}{
  an index vector for the rows of
  \code{mf = model.frame(object, "combined")} or a language object
  evaluating to such a vector.
  Only time series corresponding to indexed rows are plotted and
  only fitting windows corresponding to indexed rows are highlighted;
  the default is to plot all series and to highlight all windows.
  Evaluation of language objects follows the implementation of
  \code{\link{subset.data.frame}}.}
\item{order}{
  a permutation of \code{seq_len(nrow(mf))} or a language object
  evaluating in \code{mf} to such a vector.
  \code{order} indicates the order in which time series are plotted;
  the default indicates the original order.}
\item{xlim, ylim}{
  numeric vectors of length 2 specifying axis limits, which are
  recycled for all plots.
  If \code{time_as = "Date"}, then \code{xlim} can instead be a
  \code{\link{Date}}, \code{\link{POSIXct}}, or \code{\link{POSIXlt}}
  vector.}
\item{main, sub, xlab, ylab}{
  character or expression vectors or (\code{main}, \code{sub})
  language objects evaluating in \code{mf} to such vectors.
  These are used to generate plot (\code{main}, \code{sub}) and
  axis (\code{xlab}, \code{ylab}) labels.}
\item{\dots}{unused optional arguments.}
}
\value{
A data frame inheriting from class \code{plot.egf}.

If argument \code{cache} was supplied in the function call,
then this data frame is the result of augmenting \code{cache}
with any new computations.
}
\details{
Computation of fitted and predicted values and standard errors is
performed before any plots are created.  To avoid waste of computation
time, cached computations are returned \emph{even if} an error
is thrown during plotting.  To ensure that the cache is preserved,
assign the result of the function call to a name:
\preformatted{cache <- plot(x, \dots)}.
Caching functionality must be used with care, as mismatch between
\code{x} and \code{cache} will not be detected.  Constructions such
as \code{plot(y, cache = plot(x, \dots), \dots)}, where \code{x}
and \code{y} are different \code{egf} objects, should not be expected
to produce correct results.
}
\seealso{
The generic function \code{\link{plot}}.
}
\examples{
x <- egf_cache("egf-1.rds")
l <- list(legend = list(cex = 0.8),
          value = list(cex = 0.8, font = 2),
          ci = list(cex = 0.8))
control <- egf_control_plot(doubling = l)

op <- par(mar = c(3.5, 5, 5, 1))
plot(x,
     type = "interval",
     show_predict = 2L,
     show_doubling = 2L,
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
