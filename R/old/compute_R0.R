#' \loadmathjax
#' Compute the basic reproduction number
#'
#' Calculates basic reproduction numbers given initial exponential
#' growth rates and a binned generation interval distribution.
#'
#' @param r
#'   A non-negative numeric vector listing initial exponential
#'   growth rates. Alternatively, an `"egf"` object returned by
#'   [egf()].
#' @param breaks
#'   A numeric vector of length 2 or greater listing increasing
#'   break points in the support of the generation interval
#'   distribution, in reciprocal units of `r`
#'   (days if `inherits(r, "egf")`).
#' @param probs
#'   A numeric vector of length `length(breaks)-1`.
#'   `probs[i]` is the probability that the generation interval
#'   is between `breaks[i]` and `breaks[i+1]`.
#'   If `sum(probs) != 1`, then `probs` is replaced with
#'   `probs / sum(probs)`.
#' @inheritParams compute_tdoubling
#'
#' @details
#' Let \mjseqn{t_0 < \cdots < t_m} be the break points specified by
#' `breaks`, and let \mjseqn{p_1,\ldots,p_m} be the probabilities
#' specified by `probs`, defining \mjseqn{p_i} to be the probability
#' that the generation interval is in \mjseqn{\lbrack t_{i-1},t_i)}.
#'
#' Section 3d in \insertCite{WallLips07;textual}{epigrowthfit}
#' gives the basic reproduction number \mjseqn{\mathcal{R}USCORE0}
#' as a function of the initial exponential growth rate \mjseqn{r}:
#'
#' \mjsdeqn{\mathcal{R}USCORE0(r) = \left. r \middle/ \bigg\lbrace \sum_{i=1}^{m} \frac{p_i (e^{-r t_{i-1}} - e^{-r t_i})}{t_i - t_{i-1}} \bigg\rbrace \right.\,.}
#'
#' @return
#' The default method returns a numeric vector of length `length(r)`
#' listing basic reproduction numbers \mjseqn{\mathcal{R}USCORE0}
#' (see Details).
#'
#' The method for class `"egf"` constructs the data frame
#' `fitted(r, par = "log_r", link = FALSE, ...)`, replaces
#' variable `value` with `compute_R0(value, breaks, probs)`,
#' and returns the result.
#'
#' @examples
#' r <- seq(0, 1, by = 0.02)
#' breaks <- 0:20
#' probs <- diff(pgamma(breaks, shape = 1, scale = 2.5))
#' R0 <- compute_R0(r, breaks, probs)
#'
#' plot(r, R0, las = 1,
#'   xlab = "initial exponential growth rate, per day",
#'   ylab = "basic reproduction number"
#' )
#'
#' @references
#' \insertRef{WallLips07}{epigrowthfit}
#'
#' @seealso [compute_final_size()]
#' @export
compute_R0 <- function(r, breaks, probs, ...) {
  UseMethod("compute_R0", r)
}

#' @rdname compute_R0
#' @export
compute_R0.default <- function(r, breaks, probs, ...) {
  stop_if_not(
    is.numeric(r),
    m = "`r` must be numeric."
  )
  if (length(r) == 0L) {
    return(numeric(0L))
  }
  stop_if_not(
    is.numeric(breaks),
    length(breaks) >= 2L,
    is.finite(breaks),
    m = "`breaks` must be a finite numeric vector\nof length 2 or greater."
  )
  stop_if_not(
    diff(breaks) > 0,
    m = "`breaks` must be increasing."
  )
  stop_if_not(
    is.numeric(probs),
    length(probs) == length(breaks) - 1L,
    is.finite(probs),
    m = "`probs` must be a finite numeric vector\nof length `length(breaks)-1`."
  )
  stop_if_not(
    probs >= 0,
    any(probs > 0),
    m = "Elements of `probs` must be non-negative\nand sum to a positive number."
  )

  R0 <- rep_len(NA_real_, length(r))

  ## Degenerate cases
  if (any(r < 0, na.rm = TRUE)) {
    warning("NA returned for negative elements `r`.")
  }

  ## Limiting cases
  R0[r == 0] <- 1
  R0[r == Inf] <- Inf

  ## Usual cases
  l <- is.finite(r) & r > 0
  if (any(l)) {
    e1 <- exp(outer(breaks[-length(breaks)], -r[l]))
    e2 <- exp(outer(breaks[-1L], -r[l]))
    R0[l] <- r[l] / colSums(probs * (e1 - e2) / diff(breaks))
  }
  R0
}

#' @rdname compute_R0
#' @export
#' @importFrom stats fitted
compute_R0.egf <- function(r, breaks, probs, ...) {
  s <- c("exponential", "logistic", "richards")
  stop_if_not(
    r$curve %in% s,
    m = paste0("`r$curve` must be one of:\n", paste(dQuote(s, FALSE), collapse = ", "))
  )
  d <- fitted(r, par = "log_r", link = FALSE, ...)
  d$par <- factor(d$par, levels = "r", labels = "R0")
  d$value <- compute_R0(d$value, breaks, probs)
  d
}
