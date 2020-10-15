#' \loadmathjax
#' Compute the expected epidemic final size
#'
#' @description
#' Computes, as a function of the basic reproduction number,
#' the proportion of the population expected to be infected
#' over the course of an epidemic.
#'
#' @param R0 A numeric vector listing (non-negative) values
#'   for the basic reproduction number.
#' @param S0 A numeric vector listing values in the interval
#'   \[0,1\] for the initial susceptible proportion.
#' @param I0 A numeric vector listing values in the interval
#'   \[0,1\] for the initial infected proportion.
#'
#' @return
#' A numeric vector listing epidemic final sizes. See Details.
#'
#' The arguments are recycled up to the length of longest argument,
#' hence the length of the output is this maximum length.
#'
#' @details
#' The basic reproduction number \mjseqn{\mathcal{R}_0} defines
#' the expected epidemic final size \mjseqn{Z \in \lbrack 0,1 \rbrack}
#  implicitly through the relation
#'
#' \mjsdeqn{Z = S_0 \big(1 - e^{-\mathcal{R}_0 (Z + I_0)}\big)\,,}
#'
#' where \mjseqn{S_0,I_0 \in \lbrack 0,1 \rbrack} are the
#' proportion of the population that is susceptible and
#' infected, respectively, at the start of the epidemic.
#'
#' \insertCite{MaEarn06;textual}{epigrowthfit} discuss the history
#' and generality of this relation and show that it is valid for
#' a large class of epidemic models, including the SIR model with
#' multiple infectious stages and arbitrarily distributed stage
#' durations.
#'
#' The explicit solution for \mjseqn{Z} in terms of
#' \mjseqn{\mathcal{R}_0} involves the non-elementary
#' [Lambert W function](https://en.wikipedia.org/wiki/Lambert_W_function):
#'
#' \mjsdeqn{Z(\mathcal{R}_0) = S_0 + \frac{1}{\mathcal{R}_0} W\big\lbrack-\mathcal{R}_0 S_0 e^{-\mathcal{R}_0 (S_0 + I_0)}\big\rbrack\,.}
#'
#' @references
#' \insertRef{MaEarn06}{epigrowthfit}
#'
#' @examples
#' R0 <- seq(0, 40, by = 0.2)
#' final_size <- compute_final_size(R0)
#' plot(R0, final_size, type = "l", las = 1, ylab = "final size")
#'
#' @export
#' @importFrom emdbook lambertW
compute_final_size <- function(R0, S0 = 1, I0 = 0) {
  if (!is.numeric(R0) || length(R0) == 0) {
    stop("`R0` must be numeric and have nonzero length.")
  } else if (any(is.infinite(R0)) || isTRUE(any(R0 < 0))) {
    stop("`R0` must not contain infinite or negative values.")
  }
  if (!is.numeric(S0) || length(S0) == 0) {
    stop("`S0` must be numeric and have nonzero length.")
  } else if (isFALSE(all(S0 >= 0 & S0 <= 1))) {
    stop("Elements of `S0` must be in the interval [0,1].")
  }
  if (!is.numeric(I0) || length(I0) == 0) {
    stop("`I0` must be numeric and have nonzero length.")
  } else if (isFALSE(all(I0 >= 0 & I0 <= 1))) {
    stop("Elements of `I0` must be in the interval [0,1].")
  } else if (isFALSE(all(S0 + I0 <= 1))) {
    stop("Elements of `S0 + I0` must be in the interval [0,1].")
  }

  l <- max(length(R0), length(S0), length(I0))
  R0 <- rep(R0, length.out = l)
  S0 <- rep(S0, length.out = l)
  I0 <- rep(I0, length.out = l)
  S0 + (1 / R0) * lambertW(-R0 * S0 * exp(-R0 * (S0 + I0)))
}
