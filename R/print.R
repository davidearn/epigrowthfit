#' @export
#' @importFrom utils str capture.output
print.egf <- function(x, width = 0.9 * getOption("width"), indent = 2L, ...) {
  width <- as.integer(width)
  indent <- strrep(" ", indent)

  mean <- "f(t) - f(s)"
  if (x$model$day_of_week > 0L) {
    mean <- paste("w(s, t) *", mean)
  }
  distribution <- switch(x$model$family,
    pois = sprintf("dpois(lambda = %s)", mean),
    nbinom = sprintf("dnbinom(mu = %s, size = disp)", mean)
  )
  distribution <- paste("X(s, t) ~", distribution)
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

  get_formula_string <- function(frame, parameter, indent) {
    cl <- call("~", str2lang(parameter), formula(terms(frame))[[2L]])
    paste0(indent, capture.output(print(cl)), collapse = "\n")
  }

  cat0(wrap("This is an 'egf' object. See '?egf' for class documentation.", width = width), "\n")
  cat0("\n")

  heading("Top level nonlinear model", width = width)
  cat0("\n")
  cat0(indent, distribution, "\n")
  cat0(indent, "where\n")
  cat0(indent, cumulative, "\n")
  cat0("\n")
  str(x$model, no.list = TRUE, indent.str = indent)
  cat0("\n")

  heading("Bottom level mixed effects model", width = width)
  cat0("\n")
  w <- nchar(names(x$frame_parameters))
  formula_string <- mapply(get_formula_string,
    frame = x$frame_parameters,
    parameter = names(x$frame_parameters),
    indent = paste0(indent, strrep(" ", max(w) - w))
  )
  cat0(paste0(formula_string, "\n"))
  cat0("\n")

  heading("Data", width = width)
  dd <- data.frame(
    c(
      n <- sum(!is.na(x$frame$window)),
      N <- nlevels(x$frame$window),
      nlevels(x$frame$ts)
    ),
    c(
      sprintf("top level %s (%s)", pluralize("observation", n), pluralize("count", n)),
      sprintf("bottom level %s (fitting %s)", pluralize("observation", N), pluralize("window", N)),
      "time series"
    ),
    fix.empty.names = FALSE,
    stringsAsFactors = FALSE
  )
  cat0(paste0(indent, sub("^.", "", capture.output(print(dd, right = FALSE, quote = FALSE, row.names = FALSE, max = 6L))), "\n"))

  invisible(x)
}
