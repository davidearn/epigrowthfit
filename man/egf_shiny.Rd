\name{egf_shiny}
\alias{egf_shiny}
\title{Run Shiny applications}
\usage{
egf_shiny(app = c("toplevel", "windowselect"), ...)
}
\arguments{
\item{app}{A \link{character} string naming an application; see Details.}

\item{...}{Optional arguments to \code{\link[shiny]{runApp}}.}
}
\value{
The result of \code{\link[shiny]{runApp}}, though typically
execution is terminated before the function can return.
}
\description{
Launches a Shiny application built into \pkg{epigrowthfit}.
To halt execution, interrupt R by typing \kbd{Ctrl+C} or \kbd{Esc}.
An error is thrown if \pkg{shiny} is not installed.
}
\details{
Currently, the following applications are available:
\describe{
\item{\code{toplevel}}{
  A tool for exploring top level nonlinear models implemented
  in \pkg{epigrowthfit} and their dependence on parameters.
  Expected incidence is plotted as a function of time together
  with simulated incidence time series.
}
\item{\code{windowselect}}{
  A tool for viewing collections of time series and interactively
  selecting suitable fitting windows. Users upload time series
  data in \code{\link[=readRDS]{.rds}} format and download the
  results of their session in the same format. The result is
  a \link{list} containing arguments to be passed directly to
  \code{\link{egf}}. The element of the list specifying fitting
  windows can itself be uploaded to a session, allowing users
  to update or delete existing windows or define new ones.
}
}
}
\examples{
\dontrun{
egf_shiny("toplevel")
egf_shiny("windowselect")
}

}
