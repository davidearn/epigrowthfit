#' @export
#' @importFrom utils str
print.egf <- function(x, ...) {
  mean <- "f(t) - f(s)"
  if (x$model$day_of_week > 0L) {
    mean <- paste0("w(s, t) * ", mean)
  }
  distribution <- switch(x$model$family,
    pois = sprintf("dpois(lambda = %s)", mean),
    nbinom = sprintf("dnbinom(mu = %s, size = disp)", mean)
  )
  distribution <- paste0("X(s, t) ~ ", distribution)
  cumulative <- switch(x$model$curve,
    exponential = "c0 * exp(r * t)",
    subexponential = "c0 * (1 + alpha * (1 - p) * t / c0^(1 - p))^(1 / (1 - p))",
    gompertz = "K * exp(-exp(-alpha * (t - tinfl)))",
    logistic = "K / (1 + exp(-r * (t - tinfl)))",
    richards = "K / (1 + a * exp(-a * r * (t - tinfl)))^(1 / a)"
  )
  if (x$model$excess) {
    cumulative <- paste("b * t +", cumulative)
  }
  cumulative <- paste("f(t) =", cumulative)

  width <- as.integer(0.9 * getOption("width"))
  cat(wrap("This is an 'egf' object. See '?egf' for class details.", width = width), "\n")
  cat("\n")
  h <- "Top level nonlinear model"
  cat(h, strrep(".", max(0L, width - nchar(h))), "\n")
  cat("\n")
  cat(" ", distribution, "\n")
  cat(" ", "where\n")
  cat(" ", cumulative, "\n")
  cat("\n")
  str(x$model, no.list = TRUE, indent.str = "  ")
  cat("\n")
  h <- "Bottom level mixed effects model"
  cat(h, strrep(".", max(0L, width - nchar(h))), "\n")
  cat("\n")
  f <- function(x, s) deparse1(call("~", str2lang(s), formula(terms(x))[[2L]]))
  bottom <- mapply(f, x$frame_parameters, names(x$frame_parameters))
  w <- nchar(names(x$frame_parameters))
  cat(paste(strrep(" ", 1L + max(w) - w), bottom, "\n", collapse = ""))
  cat("\n")
  h <- "Data"
  cat(h, strrep(".", max(0L, width - nchar(h))), "\n")
  cat("\n")
  n <- sum(!is.na(x$frame$window))
  N <- nlevels(x$frame$window)
  M <- nlevels(x$frame$ts)
  w <- 1L + as.integer(log10(max(n, N, M)))
  cat(sprintf("%*d top level observation%s\n", 2L + w, n, if (n > 0L) "s" else ""))
  cat(sprintf("%*d bottom level observation%s (fitting window%s)\n", 2L + w, N, if (N > 0L) "s" else "", if (N > 0L) "s" else ""))
  cat(sprintf("%*d time series\n", 2L + w, M))
  invisible(x)
}
