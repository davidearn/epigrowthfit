#' @export
#' @importFrom utils str capture.output
print.egf <- function(x, width = 0.9 * getOption("width"), indent = 2L, ...) {
  width <- as.integer(width)
  indent <- strrep(" ", indent)

  ## Set up ===========================================================

  ## Top level nonlinear model
  top_mean <- "f(t) - f(s)"
  if (x$model$day_of_week > 0L) {
    top_mean <- paste("w(s, t) *", top_mean)
  }
  top_distribution <- switch(x$model$family,
    pois = sprintf("dpois(lambda = %s)", top_mean),
    nbinom = sprintf("dnbinom(mu = %s, size = disp)", top_mean)
  )
  top_distribution <- paste("X(s, t) ~", top_distribution)
  top_cumulative <- switch(x$model$curve,
    exponential = "c0 * exp(r * t)",
    subexponential = "c0 * (1 + alpha * (1 - p) * t / c0^(1 - p))^(1 / (1 - p))",
    gompertz = "K * exp(-exp(-alpha * (t - tinfl)))",
    logistic = "K / (1 + exp(-r * (t - tinfl)))",
    richards = "K / (1 + a * exp(-a * r * (t - tinfl)))^(1 / a)"
  )
  if (x$model$excess) {
    top_cumulative <- paste0("b * t + ", top_cumulative)
  }
  top_cumulative <- paste0("f(t) = ", top_cumulative)
  top_lines <- paste0(indent, c(top_distribution, "where", top_cumulative))

  ## Bottom level mixed effects model
  f <- function(frame, parameter, indent) {
    cl <- call("~", str2lang(parameter), formula(terms(frame))[[2L]])
    lines <- capture.output(print(cl))
    paste0(indent, lines)
  }
  nch <- nchar(names(x$frame_parameters))
  offset <- strrep(" ", max(nch) - nch)
  bottom_lines <- unlist(Map(f,
    frame = x$frame_parameters,
    parameter = names(x$frame_parameters),
    indent = paste0(indent, offset)
  ))

  ## Data
  n <- c(
    sum(!is.na(x$frame$window)),
    nlevels(x$frame$window),
    nlevels(x$frame$ts)
  )
  ltext <- format(sprintf("%d", n), justify = "right")
  rtext <- c(
    sprintf("top level %s (%s)", pluralize("observation", n[1L]), pluralize("count", n[1L])),
    sprintf("bottom level %s (fitting %s)", pluralize("observation", n[2L]), pluralize("window", n[2L])),
    "time series"
  )
  data_lines <- paste0(indent, ltext, " ", rtext)

  ## Print ============================================================

  cat(wrap("This is an 'egf' object. See '?egf' for class documentation.", width = width), "\n", sep = "")
  cat("\n")

  heading("Top level nonlinear model", width = width)
  cat("\n")
  writeLines(top_lines)
  cat("\n")
  str(x$model, no.list = TRUE, indent.str = indent, give.head = FALSE)
  cat("\n")

  heading("Bottom level mixed effects model", width = width)
  cat("\n")
  writeLines(bottom_lines)
  cat("\n")

  heading("Data", width = width)
  cat("\n")
  writeLines(data_lines)

  invisible(x)
}
