egf_shiny <-
function(app = c("toplevel", "windowselect"), ...) {
	app <- match.arg(app)
	dep <- list(toplevel = c("shiny"),
	            windowselect = c("shiny", "shinyFeedback"))
	for (package in dep[[app]]) {
		suggest(package)
	}
	dirname <- system.file("shiny", app, package = "epigrowthfit",
	                       mustWork = TRUE)
	shiny::runApp(dirname, ...)
}
