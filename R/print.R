#' @export
#' @importFrom utils str capture.output
print.egf <- function(x, width = 0.9 * getOption("width"), indent = 2L, ...) {
  width <- as.integer(width)

  ## Set up ===========================================================

  ## Top level nonlinear model
  mu <- "f(t) - f(s)"
  if (x$model$day_of_week > 0L) {
    mu <- paste("w(s, t) *", mu)
  }
  lines_top <- character(3L)
  lines_top[1L] <- switch(x$model$family,
    pois = sprintf("dpois(lambda = %s)", mu),
    nbinom = sprintf("dnbinom(mu = %s, size = disp)", mu)
  )
  lines_top[1L] <- paste("X(s, t) ~", lines_top[1L])
  lines_top[2L] <- "where"
  lines_top[3L] <- switch(x$model$curve,
    exponential = "c0 * exp(r * t)",
    subexponential = "c0 * (1 + alpha * (1 - p) * t / c0^(1 - p))^(1 / (1 - p))",
    gompertz = "K * exp(-exp(-alpha * (t - tinfl)))",
    logistic = "K / (1 + exp(-r * (t - tinfl)))",
    richards = "K / (1 + a * exp(-a * r * (t - tinfl)))^(1 / a)"
  )
  if (x$model$excess) {
    lines_top[3L] <- paste0("b * t + ", lines_top[3L])
  }
  lines_top[3L] <- paste0("f(t) = ", lines_top[3L])

  ## Bottom level mixed effects model
  f <- function(top) {
    cl <- call("~", str2lang(top), formula(x, top = top)[[2L]])
    capture.output(print(cl))
  }
  offset <- function(s) {
    n <- nchar(s)
    max(n) - n
  }
  top <- egf_get_names_top(x, link = TRUE)
  lines_bottom <- unlist1(lapply(top, f))

  ## Number of observations
  frame <- model.frame(x)
  n <- c(
    sum(!is.na(frame$x)),
    nlevels(frame$window),
    nlevels(frame$ts)
  )
  c1 <- sprintf("%d", n)
  c2 <- c(pluralize(c("observation", "fitting window"), n[1:2]), "time series")
  lines_nobs <- align(c1, c2, justify = c("r", "l"))

  ## Print ============================================================

  heading("Top level nonlinear model", width = width)
  cat("\n")
  writeLines(paste0(strrep(" ", indent), lines_top))
  cat("\n")
  str(x$model, no.list = TRUE, indent.str = strrep(" ", indent), give.head = FALSE)
  cat("\n")

  heading("Bottom level mixed effects model", width = width)
  cat("\n")
  writeLines(paste0(strrep(" ", indent + offset(top)), lines_bottom))
  cat("\n")

  heading("Data", width = width)
  cat("\n")
  writeLines(paste0(strrep(" ", indent), lines_nobs))

  invisible(x)
}

#' @export
#' @importFrom utils str capture.output
print.egf_no_fit <- print.egf
