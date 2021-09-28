#' @export
#' @importFrom stats fitted aggregate
summary.egf <- function(object, ...) {
  safe_summary <- function(x) {
    res <- unclass(summary(x))
    if (anyNA(x)) {
      return(res)
    }
    c(res, `NA's` = 0)
  }
  fo <- fitted(object)
  sfo <- aggregate(fo["estimate"], fo["top"], safe_summary)
  sfo <- `colnames<-`(t(sfo[[2L]]), as.character(sfo[[1L]]))
  if (sum(sfo[7L, ] == 0)) {
    sfo <- sfo[-7L, , drop = FALSE]
  }
  res <- list(
    convergence = object$optimizer_out$convergence,
    value = object$value,
    gradient = object$gradient,
    fitted = sfo
  )
  class(res) <- c("egf_summary", "list")
  res
}

#' @export
print.egf_summary <- function(x, width = 0.9 * getOption("width"), indent = 2L, ...) {
  width <- as.integer(width)
  indent <- strrep(" ", indent)

  heading("Negative log likelihood", width = width)
  cat("\n")
  vrag <- sprintf("%.6e", c(x$value, range(abs(x$gradient))))
  ltext <- c("convergence", "value", "range(abs(gradient))")
  rtext <- c(sprintf("%d", x$convergence), vrag[1L], paste(vrag[2L], vrag[3L]))
  writeLines(paste0(indent, format(ltext, justify = "left"), " ", rtext))
  cat("\n")

  heading("Fitted values", width = width)
  cat("\n")
  writeLines(paste0(indent, capture.output(print(x$fitted))))

  invisible(x)
}
