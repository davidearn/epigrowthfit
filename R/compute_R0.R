#' \loadmathjax
#' Compute the basic reproduction number
#'
#' @description
#' Compute the basic reproduction number corresponding
#' to an initial exponential growth rate and a binned
#' generation interval distribution.
#'
#' @param x \mjseqn{\lbrace\,r\,\rbrace}
#'   A numeric vector listing values for the initial exponential
#'   growth rate expressed per day. Alternatively, an "egf_init"
#'   or "egf" object.
#' @param breaks \mjseqn{\lbrace\,t_i\,\rbrace}
#'   A numeric vector of length 2 or greater listing increasing
#'   break points, in days, in the support of the generation
#'   interval distribution.
#' @param probs \mjseqn{\lbrace\,p_i\,\rbrace}
#'   A numeric vector with length `length(breaks)-1`. `probs[i]`
#'   is the probability that the generation interval is between
#'   `break[i]` days and `break[i+1]` days. Replaced with
#'   `probs / sum(probs)` in the event that `sum(probs) != 1`.
#'
#' @return
#' The method for class "numeric" returns a numeric vector
#' of length `length(x)`, whose `i`th element is the basic
#' reproduction number corresponding to initial exponential
#' growth rate `x[i]`. See Details.
#'
#' The method for class "egf_init" applies the method for
#' class "numeric" to `x$theta_init[["r"]]`.
#'
#' The method for class "egf" applies the method for
#' class "numeric" to `x$theta_fit[["r"]]`.
#'
#' @details
#' Let \mjseqn{t_0 < \cdots < t_m} be the break points specified
#' by `breaks`. For \mjseqn{i \in \lbrace 1,\ldots,m \rbrace},
#' let \mjseqn{p_i} be the probability specified by `probs[i]`
#' that the generation interval is in the interval
#' \mjseqn{\lbrace t_{i-1},t_i \rbrace}.
#' Section 3(d) in \insertCite{WallLips07;textual}{epigrowthfit}
#' gives the basic reproduction number \mjseqn{\mathcal{R}USCORE0}
#' as a function of the initial exponential growth rate \mjseqn{r}:
#'
#' \mjsdeqn{\mathcal{R}USCORE0(r) = \left. r \middle/ \bigg\lbrace \sum_{i=1}^{m} \frac{p_i (e^{-r t_{i-1}} - e^{-r t_i})}{t_i - t_{i-1}} \bigg\rbrace \right.\,.}
#'
#' @references
#' \insertRef{WallLips07}{epigrowthfit}
#'
#' @seealso [compute_final_size()] for computing the expected epidemic
#'   final size as a function of the basic reproduction number
#'
#' @examples
#' data(plague_latent_period)
#' lat <- plague_latent_period$relfreq
#' m <- length(lat)
#'
#' data(plague_infectious_period)
#' inf <- plague_infectious_period$relfreq
#' n <- length(inf)
#'
#' r <- seq(0, 1, by = 0.02)
#' breaks <- 1:(m+n)
#' probs <- diff(pgi(breaks, lat, inf))
#'
#' R0 <- compute_R0(r, breaks, probs)
#' plot(r, R0, las = 1,
#'   xlab = "initial exponential growth rate, per day",
#'   ylab = "basic reproduction number"
#' )
#'
#' @export
compute_R0 <- function(x, breaks, probs) {
  if (!is.numeric(breaks) || length(breaks) < 2) {
    stop("`breaks` must be numeric and have length 2 or greater.")
  } else if (!isTRUE(all(breaks >= 0))) {
    stop("Elements of `breaks` must be non-negative.")
  } else if (!all(diff(breaks) > 0)) {
    stop("`breaks` must be increasing.")
  }
  if (!is.numeric(probs) || length(probs) != length(breaks) - 1) {
    stop("`probs` must be numeric with length `length(breaks)-1`.")
  } else if (!all(is.finite(probs)) || any(probs < 0)) {
    stop("`probs` must not contain missing, infinite, or negative values.")
  } else if (all(probs == 0)) {
    stop("`probs` must have at least one positive element.")
  }

  UseMethod("compute_R0", x)
}

#' @rdname compute_R0
#' @export
compute_R0.numeric <- function(x, breaks, probs) {
  if (length(x) > 1) {
    sapply(x, compute_R0, breaks = breaks, probs = probs)
  } else {
    e1 <- exp(-x * breaks[-length(breaks)])
    e2 <- exp(-x * breaks[-1])
    d <- diff(breaks)
    x / sum(probs * (e1 - e2) / d)
  }
}

#' @rdname compute_R0
#' @export
compute_R0.egf_init <- function(x, breaks, probs) {
  x <- x$theta_init[["r"]]
  compute_R0(x, breaks, probs)
}

#' @rdname compute_R0
#' @export
compute_R0.egf <- function(x, breaks, probs) {
  x <- x$theta_fit[["r"]]
  compute_R0(x, breaks, probs)
}
