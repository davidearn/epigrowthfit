#' Compute the basic reproduction number
#'
#' Computes the basic reproduction number as a function of the
#' initial exponential growth rate, conditional on a binned
#' generation interval distribution.
#'
#' @param r
#'   A non-negative \link{numeric} vector listing initial exponential
#'   growth rates.
#' @param breaks
#'   An increasing \link{numeric} vector of length 2 or greater listing
#'   break points in the support of the generation interval distribution,
#'   in reciprocal units of \code{r}.
#' @param probs
#'   A \link{numeric} vector of length \code{\link{length}(breaks)-1}.
#'   \code{probs[i]} is the probability that the generation interval
#'   is between \code{breaks[i]} and \code{breaks[i+1]}.
#'   If \code{\link{sum}(probs) != 1}, then \code{probs} is replaced
#'   with \code{probs / \link{sum}(probs)}.
#'
#' @section Computation:
#' For an initial exponential growth rate \code{r},
#' the basic reproduction number \code{R0} is computed as
#'
#' \code{R0(r) = r / sum(probs * (exp(-r * breaks[-n]) - exp(-r * breaks[-1L])) / (breaks[-1L] - breaks[-n]))},
#'
#' where \code{n = \link{length}(breaks)}.
#'
#' @return
#' A numeric vector of length \code{\link{length}(r)} listing
#' basic reproduction numbers.
#'
#' @examples
#' r <- seq(0, 1, by = 0.02)
#' breaks <- 0:20
#' probs <- diff(pgamma(breaks, shape = 1, scale = 2.5))
#' R0 <- compute_R0(r, breaks, probs)
#'
#' plot(r, R0, las = 1,
#'   xlab = "initial exponential growth rate",
#'   ylab = "basic reproduction number"
#' )
#'
#' @references
#' Wallinga J, Lipsitch M. How generation intervals shape the relationship
#' between growth rates and reproductive numbers. Proc R Soc Lond B Biol Sci.
#' 2007;274:599--604.
#'
#' @seealso \code{\link{compute_final_size}}
#' @export
compute_R0 <- function(r, breaks, probs) {
  stop_if_not(
    is.numeric(r),
    m = "`r` must be numeric."
  )
  if (length(r) == 0L) {
    return(numeric(0L))
  }
  check_probs(probs)
  probs <- probs / sum(probs)
  stop_if_not(
    is.numeric(breaks),
    length(breaks) == length(probs) + 1L,
    is.finite(breaks),
    m = "`breaks` must be a finite numeric vector of length `length(probs)+1`."
  )

  R0 <- rep_len(NA_real_, length(r))

  ## Degenerate cases
  if (any(r < 0, na.rm = TRUE)) {
    warning("NA returned for negative elements of `r`.")
  }

  ## Limiting cases
  R0[r == 0] <- 1
  R0[r == Inf] <- Inf

  ## Usual cases
  ok <- is.finite(r) & r > 0
  if (any(ok)) {
    n <- length(breaks)
    e1 <- exp(tcrossprod(breaks[-n], -r[ok]))
    e2 <- exp(tcrossprod(breaks[-1L], -r[ok]))
    R0[ok] <- r[ok] / colSums(probs * (e1 - e2) / (breaks[-1L] - breaks[-n]))
  }
  R0
}

check_probs <- function(probs) {
  s <- deparse(substitute(probs))
  stop_if_not(
    is.numeric(probs),
    length(probs) > 0L,
    m = sprintf("`%s` must be a numeric vector of nonzero length.", s),
    n = 2L
  )
  stop_if_not(
    is.finite(probs),
    m = sprintf("`%s` must be finite.", s),
    n = 2L
  )
  stop_if_not(
    probs >= 0,
    m = sprintf("`%s` be non-negative.", s),
    n = 2L
  )
  stop_if_not(
    any(probs > 0),
    m = sprintf("`%s` must sum to a positive number.", s),
    n = 2L
  )
  invisible(NULL)
}
