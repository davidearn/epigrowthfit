#' Compute expected epidemic final size
#'
#' Computes the proportion of a population expected to be infected over
#' the course of an epidemic, as a function of the basic reproduction number.
#'
#' @param R0
#'   A \link{numeric} vector listing non-negative values for the basic
#'   reproduction number.
#' @param S0,I0
#'   \link[=numeric]{Numeric} vectors listing values in the interval [0,1]
#'   for the proportions of the population that are susceptible and infected,
#'   respectively, at the start of the epidemic. (Hence \code{S0 + I0} should
#'   not exceed 1.)
#'
#' @return
#' A \link{numeric} vector listing values in the interval [0,1]
#' for the expected epidemic final size.
#'
#' @details
#' \code{R0}, \code{S0}, and \code{I0} are recycled to length
#' \code{\link{max}(\link{lengths}(\link{list}(R0, S0, I0)))}.
#'
#' @section Computation:
#' The basic reproduction number \code{R0} defines the expected
#' epidemic final size \code{Z} through the implicit relation
#'
#' \code{Z = S0 * (1 - exp(-R0 * (Z + I0)))}.
#'
#' \code{Z} can be expressed as an explicit function of \code{R0} using the
#' \href{https://en.wikipedia.org/wiki/Lambert_W_function}{Lambert W function}:
#'
#' \code{Z(R0, S0, I0) = S0 + (1 / R0) * lambertW(-R0 * S0 * exp(-R0 * (S0 + I0)))}.
#'
#' \code{compute_final_size} evaluates this function at each supplied
#' \code{(R0, S0, I0)} triple.
#'
#' @references
#' Ma J, Earn DJD. Generality of the final size formula for an epidemic
#' of a newly invading infectious disease. Bull Math Biol. 2006;68:679–-702.
#'
#' @examples
#' R0 <- 10^seq(-3, 1, length.out = 150L)
#' final_size <- compute_final_size(R0)
#'
#' plot(R0, final_size, type = "l", las = 1,
#'   xlab = "basic reproduction number",
#'   ylab = "final size"
#' )
#'
#' @seealso \code{\link{compute_R0}}
#' @export
compute_final_size <- function(R0, S0 = 1, I0 = 0) {
  stop_if_not(
    requireNamespace("emdbook", quietly = TRUE),
    m = wrap(
      "`compute_final_size` depends on function `emdbook::lambertW`, ",
      "but package `emdbook` is not installed. Install it by running ",
      "`install.packages(\"emdbook\")`, then try again."
    )
  )
  stop_if_not(
    is.numeric(R0),
    is.numeric(S0),
    is.numeric(I0),
    m = "`R0`, `S0`, and `I0` must be numeric."
  )
  len <- lengths(list(R0, S0, I0))
  if (min(len) == 0L) {
    return(numeric(0L))
  }

  n <- max(len)
  R0 <- rep_len(R0, n)
  S0 <- rep_len(S0, n)
  I0 <- rep_len(I0, n)
  fs <- rep_len(NA_real_, n)

  ## Degenerate cases
  is_bad_triple <- (is.na(R0) | is.na(S0) | is.na(I0) |
                      R0 < 0 | S0 < 0 | I0 < 0 | S0 + I0 > 1)
  if (any(is_bad_triple)) {
    warning("NA returned for invalid (R0, S0, I0) triples.")
  }

  ## Limiting cases
  l1 <- !is_bad_triple & R0 == 0
  l2 <- !is_bad_triple & R0 == Inf
  fs[l1] <- 0
  fs[l2] <- S0[l2]
  ## Usual cases
  l3 <- !(is_bad_triple | l1 | l2)
  if (any(l3)) {
    fs[l3] <- S0[l3] + (1 / R0[l3]) * emdbook::lambertW(-R0[l3] * S0[l3] * exp(-R0[l3] * (S0[l3] + I0[l3])))
  }
  fs
}
