summary.egf <-
function(object, ...) {
	summary1 <-
	function(x) {
		ans <- unclass(summary(x))
		if (anyNA(x)) ans else c(ans, "NA's" = 0)
	}
	fo <- fitted(object)
	sfo <- aggregate(fo["estimate"], fo["top"], summary1)
	sfo <- `colnames<-`(t(sfo[[2L]]), as.character(sfo[[1L]]))
	if (max(0, sfo["NA's", ]) == 0)
		sfo <- sfo[-match("NA's", rownames(sfo), 0L), , drop = FALSE]
	ans <- list(fitted = sfo,
	            convergence = object[["optimizer_out"]][["convergence"]],
	            value = object[["value"]],
	            gradient = object[["gradient"]],
	            hessian = object[["hessian"]])
	class(ans) <- "summary.egf"
	ans
}

print.summary.egf <-
function(x, width = 0.9 * getOption("width"), indent = 2L, ...) {
	width <- as.integer(width)
	indent <- strrep(" ", indent)

	heading("Fitted values", width = width)
	cat("\n")
	writeLines(paste0(indent, capture.output(print(x[["fitted"]]))))
	cat("\n")

	heading("Negative log marginal likelihood", width = width)
	cat("\n")
	c1 <-
	    c("convergence",
	      "value",
	      "range(abs(gradient))",
	      "pos. def. Hessian")
	c2 <-
	    c(as.character(x[["convergence"]]),
	      sprintf("%.6e", x[["value"]]),
	      paste0(sprintf("%.6e", range(abs(x[["gradient"]]))), collapse = " "),
	      as.character(x[["hessian"]]))
	writeLines(paste0(indent, align(c1, c2, justify = "l")))

	invisible(x)
}
