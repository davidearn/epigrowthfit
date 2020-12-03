#' \loadmathjax
#' Compute the expected epidemic final size
#'
#' @description
#' Computes, as a function of the basic reproduction number,
#' the proportion of the population expected to be infected
#' over the course of an epidemic.
#'
#' @param R0
#'   A numeric vector listing non-negative values
#'   for the basic reproduction number.
#' @param S0
#'   A numeric vector listing values in the interval
#'   \[0,1\] for the initial susceptible proportion.
#' @param I0
#'   A numeric vector listing values in the interval
#'   \[0,1\] for the initial infected proportion.
#'
#' @return
#' A numeric vector of length `len` listing epidemic final sizes,
#' where `len` is the length of the longest argument. (Arguments
#' are recycled up to this length.) See Details.
#'
#' @details
#' The basic reproduction number \mjseqn{\mathcal{R}USCORE0} defines
#' the expected epidemic final size \mjseqn{Z \in \lbrack 0,1 \rbrack}
#  implicitly through the relation
#'
#' \mjsdeqn{Z = S_0 \big(1 - e^{-\mathcal{R}USCORE0 (Z + I_0)}\big)\,,}
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
#' \mjseqn{\mathcal{R}USCORE0} involves the non-elementary
#' [Lambert W function](https://en.wikipedia.org/wiki/Lambert_W_function):
#'
#' \mjsdeqn{Z(\mathcal{R}USCORE0) = S_0 + \frac{1}{\mathcal{R}USCORE0} W\big\lbrack-\mathcal{R}USCORE0 S_0 e^{-\mathcal{R}USCORE0 (S_0 + I_0)}\big\rbrack\,.}
#'
#' @references
#' \insertRef{MaEarn06}{epigrowthfit}
#'
#' @examples
#' R0 <- 10^seq(-3, log10(10), length.out = 150)
#' final_size <- compute_final_size(R0)
#' plot(R0, final_size, type = "l", las = 1,
#'   xlab = "basic reproduction number",
#'   ylab = "final size"
#' )
#'
#' @export
#' @importFrom emdbook lambertW
compute_final_size <- function(R0, S0 = 1, I0 = 0) {
  for (a in c("R0", "S0", "I0")) {
    x <- get(a, inherits = FALSE)
    check(x,
      what = "numeric",
      len = c(1, Inf),
      sprintf("`%s` must be numeric and have nonzero length.", a)
    )
    if (a == "R0") {
      check(x,
        val = c(0, Inf),
        sprintf("Elements of `%s` must be non-negative.", a)
      )
    } else if (a %in% c("S0", "I0")) {
      check(x,
        val = c(0, 1),
        sprintf("Elements of `%s` must be in the interval [0,1].", a)
      )
    }
  }
  len <- max(length(R0), length(S0), length(I0))
  R0 <- rep(R0, length.out = len)
  S0 <- rep(S0, length.out = len)
  I0 <- rep(I0, length.out = len)
  check(S0 + I0,
    val = c(0, 1),
    sprintf("Elements of `S0 + I0` must be in the interval [0,1].")
  )

  out <- S0 + (1 / R0) * lambertW(-R0 * S0 * exp(-R0 * (S0 + I0)))
  out <- ifelse(R0 == 0, 0, out)
  out <- ifelse(is.infinite(R0), S0, out)
  out
}
