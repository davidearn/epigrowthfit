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
  nll <- c(sprintf("%d", x$convergence), vrag[1L], paste(vrag[2L], vrag[3L]))
  dim(nll) <- c(3L, 1L)
  dimnames(nll) <- list(c("convergence", "value", "range(abs(gradient))"), "")

  heading("Negative log likelihood", width = width)
  cat0(paste0(indent, capture.output(print(nll, right = FALSE, quote = FALSE, max = 3L)), "\n"))
  cat0("\n")
  heading("Fitted values", width = width)
  cat("\n")
  cat0(paste0(indent, capture.output(print(x$fitted)), "\n"))
  invisible(x)
}
