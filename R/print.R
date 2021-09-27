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

  cat(wrap("This is an 'egf' object. See '?egf' for class documentation.", width = width), "\n", sep = "")
  cat("\n")

  heading("Top level nonlinear model", width = width)
  cat("\n")
  writeLines(paste0(indent, c(distribution, "where", cumulative)))
  cat("\n")
  str(x$model, no.list = TRUE, indent.str = indent, give.head = FALSE)
  cat("\n")

  heading("Bottom level mixed effects model", width = width)
  cat("\n")
  get_formula_lines <- function(frame, parameter, indent) {
    cl <- call("~", str2lang(parameter), formula(terms(frame))[[2L]])
    lines <- capture.output(print(cl))
    paste0(indent, lines)
  }
  w <- nchar(names(x$frame_parameters))
  writeLines(unlist(Map(get_formula_lines,
    frame = x$frame_parameters,
    parameter = names(x$frame_parameters),
    indent = paste0(indent, strrep(" ", max(w) - w))
  )))
  cat("\n")

  heading("Data", width = width)
  cat("\n")
  counts <- c(
    n <- sum(!is.na(x$frame$window)),
    N <- nlevels(x$frame$window),
    nlevels(x$frame$ts)
  )
  left <- sprintf("%d", counts)
  right <- c(
    sprintf("top level %s (%s)", pluralize("observation", n), pluralize("count", n)),
    sprintf("bottom level %s (fitting %s)", pluralize("observation", N), pluralize("window", N)),
    "time series"
  )
  w <- nchar(left)
  writeLines(paste0(indent, strrep(" ", max(w) - w), left, " ", right))

  invisible(x)
}
