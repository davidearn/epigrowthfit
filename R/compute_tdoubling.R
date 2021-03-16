#' Compute doubling time
#'
#' Calculates `log(2) / r` for non-negative `r`. This is the
#' doubling time corresponding to exponential growth rate `r`.
#'
#' @param r
#'   A non-negative numeric vector listing exponential growth rates.
#'   Alternatively, an `"egf"` object returned by [egf()].
#' @param per
#'   A positive integer indicating that `r` is a rate per `per` days,
#'   in which case the result is printed with units. Use the default
#'   (`NULL`) if `r` is unitless.
#' @param ...
#'   (For `inherits(r, "egf")` only.)
#'   Optional arguments passed to [fitted.egf()],
#'   presumably `subset` and `append`.
#'
#' @details
#' The method for class `"egf"` will issue an error if `r$curve`
#' is not an element of `c("exponential", "logistic", "richards")`.
#'
#' @return
#' The default method returns `log(2) / r` with `"tdoubling"`
#' prepending to its class. `per` is retained as an attribute
#' for use by `print.tdoubling()`.
#'
#' The method for class `"egf"` constructs the data frame
#' `fitted(r, par = "log_r", link = FALSE, ...)`, replaces
#' variable `value` with `compute_tdoubling(value, per = 1L)`,
#' and returns the result. Note that `value` always has units
#' of reciprocal days, corresponding to `per = 1L`.
#'
#' @examples
#' r <- 10^(-2:0)
#' td <- compute_tdoubling(r)
#' all.equal(td, log(2) / r)
#'
#' ## `per` attribute affects printing
#' for (i in c(1, 2, 7, 8, 365, 366)) {
#'   attr(td, "per") <- i
#'   print(td)
#' }
#'
#' @export
compute_tdoubling <- function(r, ...) {
  UseMethod("compute_tdoubling", r)
}

#' @rdname compute_tdoubling
#' @export
compute_tdoubling.default <- function(r, per = NULL, ...) {
  stop_if_not(
    is.numeric(r),
    m = "`r` must be numeric."
  )
  if (!is.null(per)) {
    stop_if_not_positive_integer(per)
  }
  if (any(r < 0, na.rm = TRUE)) {
    r[r < 0] <- NA
    warning("NA returned for negative `r`.")
  }
  td <- log(2) / r
  class(td) <- c("tdoubling", class(td))
  attr(td, "per") <- per
  td
}

#' @rdname compute_tdoubling
#' @export
#' @importFrom stats fitted
compute_tdoubling.egf <- function(r, ...) {
  s <- c("exponential", "logistic", "richards")
  stop_if_not(
    r$curve %in% s,
    m = paste0("`r$curve` must be one of:\n", paste(dQuote(s, FALSE), collapse = ", "))
  )
  d <- fitted(r, par = "log_r", link = FALSE, ...)
  d$par <- factor(d$par, levels = "r", labels = "tdoubling")
  d$value <- compute_tdoubling(d$value, per = 1L)
  d
}

#' @export
print.tdoubling <- function(x, ...) {
  per <- attr(x, "per")
  if (!is.null(per)) {
    units <- switch(as.character(per),
      `1`   = "days",
      `7`   = "weeks",
      `365` = "years (1 year = 365 days)",
      sprintf("units t = %d days", per)
    )
    cat("doubling times in ", units, ":\n\n", sep = "")
  }
  print(unclass(x))
}

