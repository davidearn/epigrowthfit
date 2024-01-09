summary.egf <-
function(object, ...) {
	safe_summary <-
	function(x) {
		res <- unclass(summary(x))
		if (anyNA(x)) {
			return(res)
		}
		c(res, `NA's` = 0)
	}
	fo <- fitted(object)
	sfo <- aggregate(fo["estimate"], fo["top"], safe_summary)
	sfo <- `colnames<-`(t(sfo[[2L]]), as.character(sfo[[1L]]))
	if (sum(sfo["NA's", ] == 0)) {
		sfo <- sfo[-match("NA's", rownames(sfo), 0L), , drop = FALSE]
	}
	res <- list(fitted = sfo,
	            convergence = object$optimizer_out$convergence,
	            value = object$value,
	            gradient = object$gradient,
	            hessian = object$hessian)
	class(res) <- "egf_summary"
	res
}

print.egf_summary <-
function(x, width = 0.9 * getOption("width"), indent = 2L, ...) {
	width <- as.integer(width)
	indent <- strrep(" ", indent)

	heading("Fitted values", width = width)
	cat("\n")
	writeLines(paste0(indent, capture.output(print(x$fitted))))
	cat("\n")

	heading("Negative log marginal likelihood", width = width)
	cat("\n")
	c1 <- c("convergence",
	        "value",
	        "range(abs(gradient))",
	        "pos. def. Hessian")
	c2 <- c(sprintf("%d", x$convergence),
	        sprintf("%.6e", x$value),
	        paste0(sprintf("%.6e", range(abs(x$gradient))), collapse = " "),
	        as.character(x$hessian))
	writeLines(paste0(indent, align(c1, c2, justify = "l")))

	invisible(x)
}
