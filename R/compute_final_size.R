#' \loadmathjax
#' Compute expected epidemic final size
#'
#' @description
#' Computes, as a function of the basic reproduction number,
#' the proportion of the population expected to be infected
#' over the course of an epidemic.
#'
#' @param R0
#'   A numeric vector listing non-negative values
#'   for the basic reproduction number.
#' @param S0,I0
#'   Numeric vectors listing values in the interval \[0,1\] for the
#'   initial susceptible and infected proportions, respectively.
#'   (Hence `S0 + I0` should not exceed 1.)
#'
#' @return
#' A numeric vector of length `max(lengths(list(R0, S0, I0)))`
#' listing epidemic final sizes (see Details).
#'
#' @details
#' The basic reproduction number \mjseqn{\mathcal{R}USCORE0} defines
#' the expected epidemic final size \mjseqn{Z \in \lbrack 0,1 \rbrack}
#  implicitly through the relation
#'
#' \mjsdeqn{Z = S_0 \big(1 - e^{-\mathcal{R}USCORE0 (Z + I_0)}\big)\,,}
#'
#' where \mjseqn{S_0,I_0 \in \lbrack 0,1 \rbrack} are the
#' proportions of the population that are susceptible and
#' infected, respectively, at the start of the epidemic.
#'
#' \insertCite{MaEarn06;textual}{epigrowthfit} discuss the history
#' and generality of this relation and show that it is valid for
#' a large class of epidemic models, including the SIR model with
#' multiple infectious stages and arbitrarily distributed stage
#' durations.
#'
#' The explicit solution for \mjseqn{Z} in terms of
#' \mjseqn{\mathcal{R}USCORE0} can be expressed using the non-elementary
#' [Lambert W function](https://en.wikipedia.org/wiki/Lambert_W_function):
#'
#' \mjsdeqn{Z(\mathcal{R}USCORE0,S_0,I_0) = S_0 + \frac{1}{\mathcal{R}USCORE0} W\big\lbrack-\mathcal{R}USCORE0 S_0 e^{-\mathcal{R}USCORE0 (S_0 + I_0)}\big\rbrack\,.}
#'
#' `R0`, `S0`, and `I0` are recycled to length
#' `max(lengths(list(R0, S0, I0)))`
#' and \mjseqn{Z} is evaluated at each resulting
#' \mjseqn{(\mathcal{R}USCORE0,S_0,I_0)} triple.
#' For invalid or incomplete triples, the value
#' is `NA`.
#'
#' @references
#' \insertRef{MaEarn06}{epigrowthfit}
#'
#' @examples
#' R0 <- 10^seq(-3, log10(10), length.out = 150L)
#' final_size <- compute_final_size(R0)
#' plot(R0, final_size, type = "l", las = 1,
#'   xlab = "basic reproduction number",
#'   ylab = "final size"
#' )
#'
#' @export
#' @importFrom emdbook lambertW
compute_final_size <- function(R0, S0 = 1, I0 = 0) {
  stop_if_not(
    is.numeric(R0),
    is.numeric(S0),
    is.numeric(I0),
    m = "`R0`, `S0`, and `I0` must be numeric."
  )
  l <- lengths(list(R0, S0, I0))
  if (min(l) == 0L) {
    return(numeric(0L))
  }

  n <- max(l)
  R0 <- rep(R0, length.out = n)
  S0 <- rep(S0, length.out = n)
  I0 <- rep(I0, length.out = n)
  is_bad_triple <- (is.na(R0) | is.na(S0) | is.na(I0) |
                      R0 < 0 | S0 < 0 | I0 < 0 | S0 + I0 > 1)
  if (any(is_bad_triple, na.rm = TRUE)) {
    warning("There were invalid or incomplete (R0,S0,I0) triples.\n",
            "Result is NA for these.")
  }

  fs <- rep.int(NA_real_, n)

  ## Limiting cases
  is_zero_R0 <- (R0 == 0)
  fs[is_zero_R0] <- 0
  is_infinite_R0 <- is.infinite(R0)
  fs[is_infinite_R0] <- S0[is_infinite_R0]

  ## Usual cases
  i <- which(!(is_bad_triple | is_zero_R0 | is_infinite_R0))
  if (length(i) > 0L) {
    fs[i] <- S0[i] + (1 / R0[i]) * lambertW(-R0[i] * S0[i] * exp(-R0[i] * (S0[i] + I0[i])))
  }
  fs
}
