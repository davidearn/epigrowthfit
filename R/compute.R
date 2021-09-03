#' Compute expected epidemic final size
#'
#' Computes the proportion of a population expected to be infected over
#' the course of an epidemic, as a function of the basic reproduction number.
#'
#' @param R0
#'   A \link{numeric} vector listing non-negative values for the basic
#'   reproduction number.
#' @param S0,I0
#'   \link[=numeric]{Numeric} vectors listing values in the interval [0,1]
#'   for the proportions of the population that are susceptible and infected,
#'   respectively, at the start of the epidemic. (Hence \code{S0 + I0} should
#'   be less than or equal to 1.)
#'
#' @return
#' A \link{numeric} vector listing values in the interval [0,1]
#' for the expected epidemic final size.
#'
#' @details
#' \code{R0}, \code{S0}, and \code{I0} are recycled to length
#' \code{\link{max}(\link{lengths}(\link{list}(R0, S0, I0)))}.
#'
#' At least one of \code{S0} and \code{I0} must be supplied.
#' If \code{S0} (\code{I0}) is supplied but not \code{I0} (\code{S0}),
#' then the latter is assigned the value of one minus the former.
#'
#' @section Computation:
#' The basic reproduction number \code{R0} defines the expected
#' epidemic final size \code{Z} through the implicit relation
#'
#' \code{Z = S0 * (1 - exp(-R0 * (Z + I0)))}.
#'
#' \code{Z} can be expressed as an explicit function of \code{R0} using the
#' \href{https://en.wikipedia.org/wiki/Lambert_W_function}{Lambert W function}:
#'
#' \code{Z(R0, S0, I0) = S0 + (1 / R0) * lambertW(-R0 * S0 * exp(-R0 * (S0 + I0)))}.
#'
#' \code{compute_final_size} evaluates this function at each supplied
#' \code{(R0, S0, I0)} triple.
#'
#' @references
#' Ma J, Earn DJD. Generality of the final size formula for an epidemic
#' of a newly invading infectious disease. Bull Math Biol. 2006;68:679–-702.
#'
#' @examples
#' R0 <- 10^seq(-3, 1, length.out = 151L)
#' final_size <- compute_final_size(R0, S0 = 1, I0 = 0)
#'
#' plot(R0, final_size, type = "l", las = 1,
#'   xlab = "basic reproduction number",
#'   ylab = "final size"
#' )
#'
#' @family epidemic parameters
#' @export
compute_final_size <- function(R0, S0, I0) {
  if (!requireNamespace("emdbook", quietly = TRUE)) {
    stop(wrap(
      "'emdbook::lambertW' is needed, but 'emdbook' is not installed. ",
      "Install it by running 'install.packages(\"emdbook\")', then try again."
    ))
  }
  stopifnot(
    is.numeric(R0),
    missing(S0) || is.numeric(S0),
    missing(I0) || is.numeric(I0)
  )
  if (missing(S0)) {
    if (missing(I0)) {
      stop("At least one of 'S0' and 'I0' must be supplied.")
    }
    S0 <- 1 - I0
  } else if (missing(I0)) {
    I0 <- 1 - S0
  }
  len <- c(length(R0), length(S0), length(I0))
  if (min(len) == 0L) {
    return(numeric(0L))
  }
  n <- max(len)
  R0 <- rep_len(R0, n)
  S0 <- rep_len(S0, n)
  I0 <- rep_len(I0, n)
  Z <- rep_len(NA_real_, n)

  ## Degenerate cases
  is_bad_triple <- (is.na(R0) | is.na(S0) | is.na(I0) |
                      R0 < 0 | S0 < 0 | I0 < 0 | S0 + I0 > 1)
  if (any(is_bad_triple)) {
    warning("NA returned for invalid (R0, S0, I0) triples.")
  }

  ## Limiting cases
  l1 <- !is_bad_triple & R0 == 0
  l2 <- !is_bad_triple & R0 == Inf
  Z[l1] <- 0
  Z[l2] <- S0[l2]

  ## Usual cases
  u <- !(is_bad_triple | l1 | l2)
  if (any(u)) {
    Z[u] <- S0[u] + emdbook::lambertW(-R0[u] * S0[u] * exp(-R0[u] * (S0[u] + I0[u]))) / R0[u]
  }
  Z
}

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
#' @family epidemic parameters
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

#' Compute doubling time
#'
#' Computes doubling times corresponding to exponential growth rates.
#'
#' @param r
#'   A non-negative \link{numeric} vector.
#' @param per
#'   A positive integer indicating that \code{r} is a rate per
#'   \code{per} days, in which case the result is printed with units.
#'   Use the default (\code{\link{NULL}}) if \code{r} is unitless.
#'
#' @return
#' \code{\link{log}(2) / r} with \code{"tdoubling"} prepended
#' to its \link{class}.
#' \code{per} is retained as an \link[=attributes]{attribute}
#' for use by \code{print.tdoubling}.
#'
#' @examples
#' r <- 10^(-2:0)
#' tdoubling <- compute_tdoubling(r)
#' all.equal(tdoubling, log(2) / r)
#'
#' ## 'per' affects printing
#' for (i in c(1, 2, 7, 8, 365, 366)) {
#'   attr(tdoubling, "per") <- i
#'   print(tdoubling)
#' }
#'
#' @family epidemic parameters
#' @export
compute_tdoubling <- function(r, per = NULL) {
  stopifnot(is.numeric(r))
  if (!is.null(per)) {
    stop_if_not_integer(per, "positive")
  }
  if (any(r < 0, na.rm = TRUE)) {
    r[r < 0] <- NA
    warning("NA returned for negative 'r'.")
  }
  tdoubling <- log(2) / r
  class(tdoubling) <- c("tdoubling", class(tdoubling))
  attr(tdoubling, "per") <- per
  tdoubling
}

#' @export
print.tdoubling <- function(x, ...) {
  per <- attr(x, "per")
  if (!is.null(per)) {
    units <- switch(as.character(per),
      `1`   = "days",
      `7`   = "weeks",
      `365` = "years (1 year = 365 days)",
      sprintf("units t = %d days", per)
    )
    cat("doubling times in ", units, ":\n\n", sep = "")
  }
  class(x) <- setdiff(class(x), "tdoubling")
  attr(x, "per") <- NULL
  NextMethod("print")
}
