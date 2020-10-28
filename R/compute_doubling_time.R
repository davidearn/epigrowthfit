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
#' The method for class "numeric" returns a "doubling_time"
#' object, a numeric vector of length `length(x)` whose
#' `i`th element is the doubling time in days corresponding
#' to initial exponential growth rate `x[i]`. See Details.
#'
#' The method for class "egf_init" applies the method for
#' class "numeric" to `x$theta_init[["r"]]`.
#'
#' The method for class "egf" applies the method for
#' class "numeric" to `x$theta_fit[["r"]]`.
#'
#' @details
#' The epidemic doubling time is the time required for cumulative
#' incidence to double at the start of an epidemic when it grows
#' roughly exponentially. The doubling time corresponding to an
#' initial exponential growth rate \mjseqn{r} is
#' \mjseqn{\frac{\log 2}{r}}.
#'
#' @examples
#' r <- seq(0, 1, by = 0.02)
#' doubling_time <- compute_doubling_time(r)
#' print(doubling_time)
#' plot(r, doubling_time, las = 1,
#'   xlab = "initial exponential growth rate, per day",
#'   ylab = "doubling time, days"
#' )
#'
#' @seealso [print.doubling_time()]
#'
#' @export
compute_doubling_time <- function(x) {
  UseMethod("compute_doubling_time", x)
}

#' @rdname compute_doubling_time
#' @export
compute_doubling_time.numeric <- function(x) {
  structure(log(2) / x, class = c("doubling_time", "numeric"))
}

#' @rdname compute_doubling_time
#' @export
compute_doubling_time.egf_init <- function(x) {
  x <- x$theta_init[["r"]]
  compute_doubling_time(x)
}

#' @rdname compute_doubling_time
#' @export
compute_doubling_time.egf <- function(x) {
  r <- x$theta_fit[["r"]]
  compute_doubling_time(r)
}
