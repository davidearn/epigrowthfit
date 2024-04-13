summary.egf <-
function(object, ...) {
	summary1 <-
	function(x) {
		ans <- unclass(summary(x))
		if (anyNA(x)) ans else c(ans, "NA's" = 0)
	}
	fo <- fitted(object)
	ans. <- aggregate(fo["estimate"], fo["top"], summary1)
	ans. <- `colnames<-`(t(ans.[[2L]]), as.character(ans.[[1L]]))
	if (max(0, ans.["NA's", ]) == 0)
		ans. <- ans.[-match("NA's", rownames(ans.), 0L), , drop = FALSE]
	ans <- list(fitted = ans.,
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
	      {
	      	r <- range(abs(x[["gradient"]]))
	      	sprintf("%.6e %.6e", r[1L], r[2L])
	      },
	      as.character(x[["hessian"]]))
	writeLines(paste0(indent, align(c1, c2, justify = "l")))

	invisible(x)
}
