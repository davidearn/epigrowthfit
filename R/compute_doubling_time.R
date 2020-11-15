#' \loadmathjax
#' Compute the epidemic doubling time
#'
#' @description
#' Compute the epidemic doubling time corresponding
#' to an initial exponential growth rate and a binned
#' generation interval distribution.
#'
#' @param r
#'   A numeric vector listing non-negative values for the initial
#'   exponential growth rate expressed per day. Alternatively, an
#'   "egf_init" or "egf" object.
#'
#' @return
#' The default method returns a numeric vector of class "doubling_time"
#' and length `length(r)` whose `i`th element is the doubling time
#' in days corresponding to initial exponential growth rate `r[i]`.
#' See Details.
#'
#' The method for class "egf_init" applies the default method to
#' \code{r = r$theta_init[["r"]]}.
#'
#' The method for class "egf" applies the default method to
#' \code{r = r$theta_fit[["r"]]}.
#'
#' @details
#' The epidemic doubling time is the time required for
#' cumulative incidence to increase by a factor of 2
#' at the start of an epidemic, when it grows roughly
#' exponentially. The doubling time corresponding to
#' an initial exponential growth rate \mjseqn{r} is
#' \mjseqn{\frac{\log 2}{r}}.
#'
#' @examples
#' r <- 10^seq(-2, 0, length.out = 150)
#' doubling_time <- compute_doubling_time(r)
#' print(doubling_time)
#' plot(r, doubling_time, las = 1,
#'   xlab = "initial exponential growth rate, per day",
#'   ylab = "doubling time, days"
#' )
#'
#' @export
compute_doubling_time <- function(r) {
  UseMethod("compute_doubling_time", r)
}

#' @rdname compute_doubling_time
#' @export
compute_doubling_time.default <- function(r) {
  check(r,
    what = "numeric",
    "`r` must be a numeric vector."
  )
  check(r,
    val = c(0, Inf),
    "Elements of `r` must be non-negative."
  )

  structure(log(2) / r, class = c("doubling_time", "numeric"))
}

#' @rdname compute_doubling_time
#' @export
compute_doubling_time.egf_init <- function(r) {
  r <- r$theta_init[["r"]]
  NextMethod("compute_doubling_time", r)
}

#' @rdname compute_doubling_time
#' @export
compute_doubling_time.egf <- function(r) {
  r <- r$theta_fit[["r"]]
  NextMethod("compute_doubling_time", r)
}
