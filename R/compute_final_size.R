#' \loadmathjax
#' Compute the expected epidemic final size
#'
#' @description
#' Computes, as a function of the basic reproduction number,
#' the proportion of the population expected to be infected
#' over the course of an epidemic in which the entire population
#' is initially susceptible.
#'
#' @param R0 A numeric vector listing (non-negative) values
#'   for the basic reproduction number.
#'
#' @return
#' A numeric vector `x` of length `length(R0)`. `x[i]` is the
#' expected epidemic final size implied by basic reproduction number
#' `R0[i]`. See Details.
#'
#' @details
#' The basic reproduction number \mjseqn{\mathcal{R}_0} defines
#' the expected epidemic final size \mjseqn{Z \in \lbrack 0,1 \rbrack}
#  implicitly through the relation
#'
#' \mjsdeqn{Z = 1 - e^{-\mathcal{R}_0 Z}\,.}
#'
#' \insertCite{MaEarn06;textual}{epigrowthfit} review and extend
#' earlier work showing that this relation is valid for a large
#' class of epidemic models, including SEIR models with arbitrarily
#' distributed infectious periods.
#'
#' The explicit solution for \mjseqn{Z} in terms of \mjseqn{\mathcal{R}_0}
#' involves the non-elementary
#' [Lambert W function](https://en.wikipedia.org/wiki/Lambert_W_function):
#'
#' \mjsdeqn{Z(\mathcal{R}_0) = 1 + \frac{1}{\mathcal{R}_0} W(-\mathcal{R}_0 e^{-\mathcal{R}_0})\,.}
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
compute_final_size <- function(R0) {
  if (missing(R0)) {
    stop("Missing argument `R0`.")
  } else if (!is.numeric(R0) || length(R0) == 0) {
    stop("`R0` must be numeric and have nonzero length.")
  } else if (any(is.infinite(R0)) || isTRUE(any(R0 < 0))) {
    stop("`R0` must not contain infinite or negative values.")
  }
  1 + (1 / R0) * lambertW(-R0 * exp(R0))
}
