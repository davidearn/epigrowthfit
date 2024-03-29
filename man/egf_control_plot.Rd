\name{egf_control_plot}
\alias{egf_control_plot}
\title{Define plot options}
\usage{
egf_control_plot(
  window,
  data,
  predict,
  asymptote,
  box,
  axis,
  title,
  tdoubling
)
}
\arguments{
\item{window}{A named list of arguments to \code{\link{rect}}
affecting the appearance of fitting windows.}

\item{data}{A named list of the form \code{list(main, short, long)}.
\code{main} is a named list of arguments to \code{\link{points}}
affecting the appearance of observed data.
\code{short} and \code{long} are alternatives to \code{main}
used for counts over intervals shorter or longer than \code{dt}
when \code{type = "interval"}.
\code{short} and \code{long} default to \code{main} (elementwise).}

\item{predict}{A named list of the form \code{list(estimate, ci)}.
\code{estimate} and \code{ci} are named lists of arguments
to \code{\link{lines}} and \code{\link{polygon}},
affecting the appearance of predicted curves and corresponding
confidence bands.
\code{ci$col} defaults to \code{estimate$col} with added transparency.}

\item{asymptote}{A named list of arguments to \code{\link{segments}}
affecting appearance of line segments drawn at
\code{y = <initial exponential growth rate>}
when \code{type = "rt"} and
\code{x$model$curve = "logistic"} or \code{"richards"}.}

\item{box}{A named list of arguments to \code{\link{box}}
affecting the appearance of the box drawn around the plot region
on the device.}

\item{axis}{A named list of the form \code{\link{list}(x, y)}.
\code{x} and \code{y} are named lists of arguments to \code{\link{axis}}
affecting the appearance of the bottom and left axes.
When \code{time_as = "Date"}, there are minor and major bottom axes.
In this case, the major axis uses a modified version of \code{x}
that tries to ensure that it is displayed below the minor axis in
a slightly larger font.}

\item{title}{A named list of the form \code{list(main, sub, xlab, ylab)}.
The elements are named lists of arguments to \code{\link{title}}
affecting the appearance of plot (sub)titles and axis labels.
\code{sub$adj} defaults to \code{main$adj}.}

\item{tdoubling}{A named list of the form \code{list(legend, estimate, ci)}.
The elements are named lists of arguments to \code{\link{mtext}}
affecting the appearance of initial doubling times printed in the
top margin.}
}
\value{
A named list.
}
\description{
Sets parameters controlling the graphical output of \code{\link{plot.egf}}.
Supplied values override package defaults
(retrievable as \code{defaults <- egf_control_plot()}),
which in turn override global defaults set via \code{\link{par}}.\cr\cr
Here, \code{x}, \code{type}, \code{time_as}, and \code{dt}
refer to the so-named arguments of \code{plot.egf}.
}
\details{
Setting an argument (or an element thereof, in the case of nested lists)
to \code{\link{NULL}} has the effect of suppressing the corresponding
plot element.
}
