#' Compute doubling time
#'
#' Calculates `log(2) / r` for non-negative `r`. This is the
#' doubling time corresponding to exponential growth rate `r`.
#'
#' @param r
#'   A numeric vector with non-negative elements listing exponential
#'   growth rates. Alternatively, an `"egf"` object returned by [egf()].
#' @param per
#'   A positive integer indicating that `r` is a rate per `per` days,
#'   in which case the result is printed with units. Use the default
#'   (`NULL`) if `r` is unitless.
#' @param subset
#'   A named list passed to [fitted.egf()].
#' @param ...
#'   Arguments passed to methods.
#'
#' @details
#' The method for class `"egf"` will issue an error if `r$curve`
#' is not an element of `c("exponential", "logistic", "richards")`.
#'
#' @return
#' The default method returns `log(2) / r`.
#'
#' The method for class `"egf"` constructs the
#' data frame `d = fitted(r, subset, link = TRUE)`,
#' appends `tdoubling = log(2) / exp(d$log_r)`,
#' and returns the result omitting extraneous variables.
#' (Here, `exp(d$log_r)` always has units of reciprocal
#' days, corresponding to `per = 1L`.)
#'
#' Unless `per = NULL`, the result has `"tdoubling"`
#' prepended to its class and `per` stored as an
#' attribute for use by `print.tdoubling()`.
#'
#' @examples
#' r <- 10^(-2:0)
#' td <- compute_tdoubling(r)
#' all.equal(td, log(2) / r)
#'
#' ## Printing depends on `per` attribute
#' class(td) <- c("tdoubling", class(td))
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
  if (any(r < 0, na.rm = TRUE)) {
    r[r < 0] <- NA
    warning("Negative elements of `r` replaced with NA.")
  }

  td <- log(2) / r
  if (!is.null(per)) {
    stop_if_not_positive_integer(per)
    attr(td, "per") <- per
    class(td) <- c("tdoubling", class(td))
  }
  td
}

#' @rdname compute_tdoubling
#' @export
#' @importFrom stats fitted
compute_tdoubling.egf <- function(r, subset, ...) {
  s <- c("exponential", "logistic", "richards")
  stop_if_not(
    r$curve %in% s,
    m = paste0(
      "`r$curve` must be one of:\n",
      paste(sprintf("\"%s\"", s), collapse = ", ")
    )
  )

  d <- fitted(r, subset, link = TRUE)
  n <- length(r$frame) - 2L # first `n` variables in `d` are factors
  out <- cbind(d[seq_len(n)], tdoubling = log(2) / exp(d$log_r))
  attr(out, "per") <- 1L
  class(out) <- c("tdoubling", class(out))
  out
}

#' @export
print.tdoubling <- function(x, ...) {
  per <- attr(x, "per")
  units <- switch(as.character(per),
    "1"   = "days",
    "7"   = "weeks",
    "365" = "years (1 year = 365 days)",
    sprintf("units t = %d days", per)
  )
  attr(x, "per") <- NULL
  class(x) <- setdiff(class(x), "tdoubling")
  cat("doubling times in ", units, ":\n\n", sep = "")
  print(x)
}

