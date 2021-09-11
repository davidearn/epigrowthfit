#' @export
#' @importFrom stats fitted aggregate
summary.egf <- function(object, ...) {
  fo <- fitted(object)
  sfo <- aggregate(fo["estimate"], fo["top"], summary)
  res <- list(
    convergence = object$optimizer_out$convergence,
    value = object$value,
    gradient = object$gradient,
    fitted = `colnames<-`(t(sfo[[2L]]), as.character(sfo[[1L]]))
  )
  class(res) <- c("egf_summary", "list")
  res
}

#' @export
print.egf_summary <- function(x, width = 0.9 * getOption("width"), indent = 2L, ...) {
  width <- as.integer(width)
  indent <- strrep(" ", indent)

  vrag <- sprintf("%.6e", c(x$value, range(abs(x$gradient))))
  dd <- data.frame(
    c(
      "convergence",
      "value",
      "range(abs(gradient))"
    ),
    c(
      sprintf("%d", x$convergence),
      vrag[1L],
      paste(vrag[2L], vrag[3L])
    ),
    fix.empty.names = FALSE,
    stringsAsFactors = FALSE
  )

  heading("Negative log likelihood", width = width)
  lines <- capture.output(print(dd, right = FALSE, quote = FALSE, row.names = FALSE, max = 6L))
  lines <- paste0(indent, sub("^.", "", lines))
  writeLines(lines)
  cat("\n")

  heading("Fitted values", width = width)
  cat("\n")
  lines <- capture.output(print(x$fitted))
  lines <- paste0(indent, lines)
  writeLines(lines)

  invisible(x)
}
