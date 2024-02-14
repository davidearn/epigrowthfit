egf_shiny <-
function(app = c("toplevel", "windowselect"), ...) {
	app <- match.arg(app)
	stopifnot(requireNamespace("shiny"))
	if (app == "windowselect")
	stopifnot(requireNamespace("shinyFeedback"))
	appdir <- system.file("shiny", app, package = "epigrowthfit",
	                      mustWork = TRUE)
	shiny::runApp(appdir, ...)
}
