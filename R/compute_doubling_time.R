#' \loadmathjax
#' Compute the epidemic doubling time
#'
#' @description
#' Compute the epidemic doubling time corresponding
#' to an initial exponential growth rate and a binned
#' generation interval distribution.
#'
#' @param x \mjseqn{\lbrace\,r\,\rbrace}
#'   A numeric vector listing values for the initial exponential
#'   growth rate expressed per day. Alternatively, an "egf_init"
#'   or "egf" object.
#'
#' @return
#' The method for class "numeric" returns a numeric vector
#' of length `length(x)`, whose `i`th element is the doubling
#' time in days corresponding to initial exponential growth
#' rate `x[i]`. See Details.
#'
#' The method for class "egf_init" applies the method for
#' class "numeric" to `x$theta0[["r"]]`.
#'
#' The method for class "egf" applies the method for
#' class "numeric" to `x$theta_hat[["r"]]`.
#'
#' @details
#' The epidemic doubling time is the time required for
#' cumulative incidence to double at the start of epidemic
#' when cumulative incidence grows roughly exponentially.
#' The doubling time corresponding to an initial exponential
#' growth rate \mjseqn{r} is \mjseqn{\frac{\log 2}{r}}.
#'
#' @examples
#' r <- seq(0, 1, by = 0.02)
#' doubling_time <- compute_doubling_time(r)
#' plot(r, doubling_time, las = 1,
#'   xlab = "initial exponential growth rate, per day",
#'   ylab = "doubling time, days"
#' )
#'
#' @export
compute_doubling_time <- function(x) {
  UseMethod("compute_doubling_time", x)
}

#' @rdname compute_doubling_time
#' @export
compute_doubling_time.numeric <- function(x) {
  if (length(x) > 1) {
    sapply(x, compute_doubling_time)
  } else {
    log(2) / x
  }
}

#' @rdname compute_doubling_time
#' @export
compute_doubling_time.egf_init <- function(x) {
  x <- x$theta0[["r"]]
  compute_doubling_time(x)
}

#' @rdname compute_doubling_time
#' @export
compute_doubling_time.egf <- function(x) {
  r <- x$theta_hat[["r"]]
  compute_doubling_time(r)
}
