#' \loadmathjax
#' Compute the basic reproduction number
#'
#' Calculates basic reproduction numbers given initial exponential
#' growth rates and a binned generation interval distribution.
#'
#' @param r
#'   A numeric vector with non-negative elements listing initial
#'   exponential growth rates \mjseqn{r}. Alternatively, an `"egf"`
#'   object returned by [egf()].
#' @param breaks
#'   A numeric vector of length 2 or greater listing increasing
#'   break points \mjseqn{t_i} in the support of the generation
#'   interval distribution, in reciprocal units of `r` (days if
#'   `r` is an `"egf"` object).
#' @param probs
#'   A numeric vector of length `length(breaks)-1`. `probs[i]`
#'   is the probability \mjseqn{p_i} that the generation interval
#'   is between `breaks[i]` and `breaks[i+1]`. If `sum(probs) != 1`,
#'   then `probs` is replaced with `probs / sum(probs)`.
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
#' `d = fitted(r, subset, link = TRUE)`, appends
#' `R0 = compute_R0.default(r = exp(d$log_r), breaks, probs)`,
#' and returns the result omitting extraneous variables.
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
#' r <- seq(0, 1, by = 0.02) # per day
#' breaks <- seq(1, m + n, by = 1) # days
#' probs <- diff(pgi(breaks, lat, inf))
#'
#' R0 <- compute_R0(r, breaks, probs)
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
    m = paste0(
      "`breaks` must be a finite numeric vector\n",
      "of length 2 or greater."
    )
  )
  stop_if_not(
    diff(breaks) > 0,
    m = "`breaks` must be increasing."
  )
  stop_if_not(
    is.numeric(probs),
    length(probs) == length(breaks) - 1L,
    is.finite(breaks),
    m = paste0(
      "`probs` must be a finite numeric vector\n",
      "of length `length(breaks)-1`."
    )
  )
  stop_if_not(
    probs >= 0,
    any(probs > 0),
    m = paste0(
      "Elements of `probs` must be non-negative\n",
      "and sum to a positive number."
    )
  )
  if (any(r < 0, na.rm = TRUE)) {
    r[r < 0] <- NA
    warning("Negative elements of `r` replaced with NA.")
  }

  R0 <- rep.int(NA_real_, length(r))

  ## Limiting cases
  is_zero_r <- (r == 0)
  R0[is_zero_r] <- 1
  is_infinite_r <- is.infinite(r)
  R0[is_infinite_r] <- Inf

  ## Usual cases
  i <- which(!(is.na(r) | is_zero_r | is_infinite_r))
  if (length(i) > 0L) {
    e1 <- exp(outer(breaks[-length(breaks)], -r[i]))
    e2 <- exp(outer(breaks[-1L], -r[i]))
    R0[i] <- r[i] / colSums(probs * (e1 - e2) / diff(breaks))
  }
  R0
}

#' @rdname compute_R0
#' @export
#' @importFrom stats fitted
compute_R0.egf <- function(r, breaks, probs, subset, ...) {
  s <- c("exponential", "logistic", "richards")
  stop_if_not(
    r$curve %in% s,
    m = paste0(
      "`r$curve` must be one of:\n",
      paste(sprintf("\"%s\"", s), collapse = ", ")
    )
  )

  d <- fitted(r, subset = subset, link = TRUE)
  cbind(
    d[vapply(d, is.factor, FALSE)],
    R0 = compute_R0(exp(d$log_r), breaks, probs)
  )
}
