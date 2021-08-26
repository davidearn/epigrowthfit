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
  stopifnot(is.numeric(r))
  if (length(r) == 0L) {
    return(numeric(0L))
  }
  stopifnot(
    is.numeric(breaks),
    length(breaks) >= 2L,
    is.finite(breaks),
    diff(breaks) > 0
  )
  stopifnot(
    is.numeric(probs),
    length(probs) == length(breaks) - 1L,
    is.finite(probs),
    probs >= 0,
    any(probs > 0)
  )
  probs <- probs / sum(probs)
  R0 <- rep_len(NA_real_, length(r))

  ## Degenerate cases
  if (any(r < 0, na.rm = TRUE)) {
    warning("NA returned for negative elements of 'r'.")
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
