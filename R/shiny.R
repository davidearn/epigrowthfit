#' Run Shiny applications
#'
#' Launches a Shiny application built into \pkg{epigrowthfit}.
#' To halt execution, interrupt R by typing \kbd{Ctrl+C} or \kbd{Esc}.
#' An error is thrown if \pkg{shiny} is not installed.
#'
#' @param app
#'   A \link{character} string naming an application; see Details.
#' @param ...
#'   Optional arguments to \code{\link[shiny]{runApp}}.
#'
#' @details
#' Currently, the following applications are available:
#' \describe{
#' \item{\code{toplevel}}{
#'   A tool for exploring top level nonlinear models implemented
#'   in \pkg{epigrowthfit} and their dependence on parameters.
#'   Expected incidence is plotted as a function of time together
#'   with simulated incidence time series.
#' }
#' \item{\code{windowselect}}{
#'   A tool for viewing collections of time series and interactively
#'   selecting suitable fitting windows. Users upload time series
#'   data in \code{\link[=readRDS]{.rds}} format and download the
#'   results of their session in the same format. The result is
#'   a \link{list} containing arguments to be passed directly to
#'   \code{\link{egf}}. The element of the list specifying fitting
#'   windows can itself be uploaded to a session, allowing users
#'   to update or delete existing windows or define new ones.
#' }
#' }
#'
#' @return
#' The result of \code{\link[shiny]{runApp}}, though typically
#' execution is terminated before the function can return.
#'
#' @examples
#' \dontrun{
#' egf_shiny("toplevel")
#' egf_shiny("windowselect")
#' }
#'
#' @export
egf_shiny <- function(app = c("toplevel", "windowselect"), ...) {
  app <- match.arg(app)
  suggests <- list(
    toplevel = c("shiny"),
    windowselect = c("shiny", "shinyFeedback")
  )
  for (s in suggests[[app]]) {
    if (!requireNamespace(s, quietly = TRUE)) {
      stop(wrap(
        "Application ", sQuote(app), " depends on uninstalled package ", sQuote(s), ". ",
        "Install it by running 'install.packages(", dQuote(s), ")', then try again."
      ))
    }
  }
  dir <- system.file("shiny", app, package = "epigrowthfit", mustWork = TRUE)
  shiny::runApp(dir, ...)
}
