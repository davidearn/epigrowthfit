\name{egf_shiny}
\alias{egf_shiny}
\title{Installed Shiny Applications}
\description{
Launches one of the \href{https://shiny.posit.co/}{Shiny} applications
installed with package \pkg{epigrowthfit}.  To halt execution, interrupt
\R{} by typing \kbd{Ctrl+C} or \kbd{Esc}.  An error is thrown if package
\pkg{shiny} is not installed.
}
\usage{
egf_shiny(app = c("toplevel", "windowselect"), \dots)
}
\arguments{
\item{app}{
  a character string naming an application.}
\item{\dots}{
  optional arguments to \code{\link[shiny]{runApp}}.}
}
\value{
The value of the \code{\link[shiny]{runApp}} call, though typically
execution is terminated before the function can return.
}
\details{
Currently, two applications are available:
\describe{
\item{\code{"toplevel"}}{
  a tool for the exploring top level nonlinear models implemented
  in \pkg{epigrowthfit} and their dependence on parameter values.
  Expected incidence is plotted as a function of time together with
  simulated incidence time series.}
\item{\code{"windowselect"}}{
  a tool for viewing collections of incidence time series and
  interactively selecting suitable fitting windows.  Users upload
  time series data in \code{\link[=saveRDS]{.rds}} format and
  download the results of their session in the same format.
  The result is a list containing arguments to be passed directly
  to \code{\link{egf}}.  The element of the list specifying
  fitting windows can itself be uploaded to a session, allowing
  users to modify or delete existing windows and define new ones.}
}
}
\examples{
\dontrun{
egf_shiny("toplevel")
egf_shiny("windowselect")
}
}
