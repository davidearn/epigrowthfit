#' @export
#' @importFrom utils str capture.output
print.egf <- function(x, width = 0.9 * getOption("width"), exdent = 2L, ...) {
  width <- as.integer(width)
  exdent <- strrep(" ", exdent)

  ## Sort out details of top level nonlinear model
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

  ## Utilities
  cat0 <- function(...) {
    cat(..., sep = "")
  }
  heading <- function(s) {
    cat0(s, " ", strrep(".", max(0L, width - 1L - nchar(s))), "\n")
  }
  pluralize <- function(s, n) {
    if (n > 1L) paste0(s, "s") else s
  }
  get_formula_string <- function(frame, parameter, indent) {
    cl <- call("~", str2lang(parameter), formula(terms(frame))[[2L]])
    paste0(indent, capture.output(print(cl)), collapse = "\n")
  }

  cat0(wrap("This is an 'egf' object. See '?egf' for class documentation.", width = width), "\n")
  cat0("\n")

  heading("Top level nonlinear model")
  cat0("\n")
  cat0(exdent, distribution, "\n")
  cat0(exdent, "where\n")
  cat0(exdent, cumulative, "\n")
  cat0("\n")
  str(x$model, no.list = TRUE, indent.str = exdent)
  cat0("\n")

  heading("Bottom level mixed effects model")
  cat0("\n")
  w <- nchar(names(x$frame_parameters))
  formula_string <- mapply(get_formula_string,
    frame = x$frame_parameters,
    parameter = names(x$frame_parameters),
    indent = paste0(exdent, strrep(" ", max(w) - w))
  )
  cat0(paste0(formula_string, "\n", collapse = ""))
  cat0("\n")

  heading("Data")
  cat0("\n")
  n <- sum(!is.na(x$frame$window))
  N <- nlevels(x$frame$window)
  M <- nlevels(x$frame$ts)
  w <- 1L + as.integer(log10(max(n, N, M)))
  cat0(exdent, sprintf("%*d top level %s\n", w, n, pluralize("observation", n)))
  cat0(exdent, sprintf("%*d bottom level %s (fitting %s)\n", w, N, pluralize("observation", N), pluralize("window", N)))
  cat0(exdent, sprintf("%*d time series\n", w, M))

  invisible(x)
}
