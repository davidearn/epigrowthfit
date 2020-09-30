#' \loadmathjax
#' Compute the basic reproduction number
#'
#' @description
#' Given an initial growth rate and a generation interval
#' distribution, compute the basic reproduction number.
#'
#' @param r A numeric vector listing values for the initial
#'   growth rate \mjseqn{r} expressed per year.
#' @param breaks A numeric vector of length 2 or greater
#'   listing increasing break points, in days, in the support
#'   of the generation interval density.
#' @param probs A numeric vector with length `length(breaks)-1`.
#'   `probs[i]` is the probability that the generation interval
#'   is between `break[i]` days and `break[i+1]` days. Replaced
#'   with `probs / sum(probs)` in the event that `sum(probs) != 1`.
#'
#' @return
#' A numeric vector of length `length(r)` listing values
#' for the basic reproduction number \mjseqn{\mathcal{R}_0}.
#' The `i`th element is conditional on initial growth rate
#' `r[i]`. See Details.
#'
#' @details
#' Let \mjseqn{a_1 < \cdots < a_\ell} be the break points
#' specified by `breaks`.
#' For \mjseqn{i \in \lbrace 1,\ldots,\ell-1 \rbrace}, let
#' \mjseqn{z_i} be the probability specified by `probs[i]`
#' that the generation interval \mjseqn{Z} is in the interval
#' \mjseqn{(a_i,a_{i+1}\rbrack}.
#'
#' Section 3(d) in \insertCite{WallLips07;textual}{epigrowthfit}
#' states that the basic reproduction number is given by
#'
#' \mjsdeqn{\out{\mathcal{R}_0 = \left. r \middle/ \bigg\lbrace \sum_{i=1}^{\ell-1} \frac{z_i (e^{-r a_i} - e^{-r a_{i+1}})}{a_{i+1} - a_i} \bigg\rbrace \right.\,,}}
#'
#' where \mjseqn{r} is the initial growth rate specified by `r`.
#'
#' @references
#' \insertRef{WallLips07}{epigrowthfit}
#'
#' @seealso [dgi()], [pgi()] for evaluating generation interval
#'   density and distribution functions assuming discrete latent
#'   and infectious period distributions
#'
#' @examples
#' data(plague_latent_period)
#' latent <- plague_latent_period$relfreq
#' m <- length(latent)
#'
#' data(plague_infectious_period)
#' infectious <- plague_latent_period$relfreq
#' n <- length(infectious)
#'
#' r <- seq(0, 1, by = 0.02)
#' breaks <- 1:(m+n)
#' probs <- diff(pgi(breaks, latent, infectious))
#'
#' R0 <- compute_R0(r, breaks, probs)
#' plot(r, R0, las = 1,
#'   xlab = "initial growth rate, per day",
#'   ylab = "basic reproduction number"
#' )
#'
#' @export
compute_R0 <- function(r, breaks, probs) {
  if (missing(r)) {
    stop("Missing argument `r`.")
  } else if (!is.numeric(r) || length(r) == 0) {
    stop("`r` must be numeric and have nonzero length.")
  }
  if (missing(breaks)) {
    stop("Missing argument `breaks`.")
  } else if (!is.numeric(breaks) || length(breaks) < 2) {
    stop("`breaks` must be numeric and have length 2 or greater.")
  } else if (!isTRUE(all(breaks >= 0))) {
    stop("Elements of `breaks` must be non-negative.")
  } else if (!all(diff(breaks) > 0)) {
    stop("`breaks` must be increasing.")
  }
  if (missing(probs)) {
    stop("Missing argument `probs`.")
  } else if (!is.numeric(probs) || length(probs) != length(breaks) - 1) {
    stop("`probs` must be numeric with length `length(breaks)-1`.")
  } else if (!all(is.finite(probs)) || any(probs < 0)) {
    stop("`probs` must not contain missing, infinite, or negative values.")
  } else if (all(probs == 0)) {
    stop("`probs` must have at least one positive element.")
  }

  if (length(r) == 1) {
    if (!is.finite(r)) {
      return(NA)
    }
    x1 <- exp(-r * breaks[-length(breaks)])
    x2 <- exp(-r * breaks[-1])
    d <- diff(breaks)
    r / sum(probs * (x1 - x2) / d)
  } else {
    sapply(r, compute_R0, breaks = breaks, probs = probs)
  }
}
