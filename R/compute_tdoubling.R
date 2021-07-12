#' Compute doubling time
#'
#' Computes \code{\link{log}(2) / r} for non-negative \code{r}, which
#' is the doubling time corresponding to exponential growth rate \code{r}.
#'
#' @param r
#'   A non-negative \link{numeric} vector listing exponential growth rates.
#' @param per
#'   A positive integer indicating that \code{r} is a rate per
#'   \code{per} days, in which case the result is printed with units.
#'   Use the default (\code{\link{NULL}}) if \code{r} is unitless.
#'
#' @return
#' \code{\link{log}(2) / r} with \code{"tdoubling"} prepended to its
#' \link{class}. \code{per} is retained as an \link[=attributes]{attribute}
#' for use by \code{print.tdoubling}.
#'
#' @examples
#' r <- 10^(-2:0)
#' tdoubling <- compute_tdoubling(r)
#' all.equal(tdoubling, log(2) / r)
#'
#' ## Attribute `per` affects printing
#' for (i in c(1, 2, 7, 8, 365, 366)) {
#'   attr(tdoubling, "per") <- i
#'   print(tdoubling)
#' }
#'
#' @export
compute_tdoubling <- function(r, per = NULL) {
  stop_if_not(
    is.numeric(r),
    m = "`r` must be numeric."
  )
  if (!is.null(per)) {
    stop_if_not_integer(per, "positive")
  }
  if (any(r < 0, na.rm = TRUE)) {
    r[r < 0] <- NA
    warning("NA returned for negative `r`.")
  }
  tdoubling <- log(2) / r
  class(tdoubling) <- c("tdoubling", class(tdoubling))
  attr(tdoubling, "per") <- per
  tdoubling
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
  class(x) <- setdiff(class(x), "tdoubling")
  attr(x, "per") <- NULL
  NextMethod("print")
}
