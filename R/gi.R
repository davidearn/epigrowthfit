#' Generation interval distribution
#'
#' Generation interval
#' density function (\code{dgi}), distribution function (\code{pgi}),
#' quantile function (\code{qgi}), and sampling (\code{rgi}).
#' Results are conditional on supplied latent and infectious period
#' distributions. It is assumed
#' \itemize{
#' \item that the latent period and infectious waiting time are independent,
#' \item that infectiousness is constant over the infectious period, and
#' \item that the latent and infectious periods are positive and integer-valued
#' (in arbitrary but like units of time).
#' }
#'
#' @param x,q
#'   A \link{numeric} vector listing generation intervals.
#' @param p
#'   A \link{numeric} vector listing probabilities.
#' @param n
#'   A non-negative integer indicating a desired sample size.
#' @param latent,infectious
#'   \link[=numeric]{Numeric} vectors such that \code{latent[i]} and
#'   \code{infectious[i]} are the probabilities that the latent and
#'   infectious periods, respectively, are \code{i} units of time.
#'   It is sufficient to supply probability weights, as both vectors
#'   are divided by their sums internally.
#'
#' @return
#' A \link{numeric} vector with length equal to the that of the
#' first argument, or length \code{n} in the case of \code{rgi}.
#'
#' @references
#' Svensson, Å. A note on generation times in epidemic models.
#' Math Biosci. 2007;208:300--11.
#'
#' @examples
#' data(plague_latent_period)
#' latent <- plague_latent_period$relfreq
#' m <- length(latent)
#'
#' data(plague_infectious_period)
#' infectious <- plague_infectious_period$relfreq
#' n <- length(infectious)
#'
#' ## Histogram of samples
#' y <- rgi(1e06, latent, infectious)
#' hist(y, breaks = seq(0, m + n + 1), freq = FALSE, las = 1,
#'   ylab = "relative frequency",
#'   main = ""
#' )
#'
#' ## Density and distribution functions
#' x <- seq(0, m + n + 1, by = 0.02)
#' fx <- dgi(x, latent, infectious)
#' Fx <- pgi(x, latent, infectious)
#' plot(x, fx, type = "l", las = 1, # consistent with histogram
#'   xlab = "generation interval",
#'   ylab = "density function"
#' )
#' plot(x, Fx, type = "l", las = 1,
#'   xlab = "generation interval",
#'   ylab = "distribution function"
#' )
#'
#' ## Quantile function
#' p <- seq(0, 1, by = 0.001)
#' qp <- qgi(p, latent, infectious)
#' plot(p, qp, type = "l", las = 1,
#'   xlab = "probability",
#'   ylab = "quantile function"
#' )
#'
#' @name gi
NULL

#' @rdname gi
#' @export
dgi <- function(x, latent, infectious) {
  stop_if_not(
    is.numeric(x),
    m = "`x` must be numeric."
  )
  if (length(x) == 0L) {
    return(x)
  }
  check_probs(latent)
  check_probs(infectious)
  latent <- latent / sum(latent)
  infectious <- infectious / sum(infectious)
  m <- length(latent)
  n <- length(infectious)

  d <- x
  d[] <- NA
  d[x < 1 | x >= m + n] <- 0

  l <- !is.na(x) & x >= 1 & x < m + n
  if (any(l)) {
    ## Take advantage of the fact that the generation interval density
    ## is constant on intervals [i,i+1)
    x_floor <- floor(x[l])
    x_floor_unique <- unique(x_floor)
    k <- match(x_floor, x_floor_unique)

    ij_dim <- c(m, length(x_floor_unique))
    i <- .row(ij_dim)
    j <- .col(ij_dim)

    d[l] <- colSums(latent * diwt(x_floor_unique[j] - i, infectious = infectious))[k]
  }
  d
}

#' @rdname gi
#' @export
#' @importFrom stats approx
pgi <- function(q, latent, infectious) {
  stop_if_not(
    is.numeric(q),
    m = "`q` must be numeric."
  )
  if (length(q) == 0L) {
    return(q)
  }
  m <- length(latent)
  n <- length(infectious)

  p <- q
  p[] <- NA
  p[q <= 1] <- 0
  p[q >= m + n] <- 1
  l <- !is.na(q) & q > 1 & q < m + n
  if (any(l)) {
    p[l] <- approx(
      x = seq_len(m + n),
      y = c(0, cumsum(dgi(seq_len(m + n - 1L), latent = latent, infectious = infectious))),
      xout = q[l]
    )$y
  }
  p
}

#' @rdname gi
#' @export
#' @importFrom stats approx
qgi <- function(p, latent, infectious) {
  stop_if_not(
    is.numeric(p),
    m = "`p` must be numeric."
  )
  if (length(p) == 0L) {
    return(p)
  }
  m <- length(latent)
  n <- length(infectious)

  q <- p
  q[] <- NA
  q[p < 0 | p > 1] <- NaN
  q[p == 0] <- -Inf
  q[p == 1] <- m + n
  l <- !is.na(p) & p > 0 & p < 1
  if (any(l)) {
    q[l] <- approx(
      x = c(0, cumsum(dgi(seq_len(m + n - 1L), latent = latent, infectious = infectious))),
      y = seq_len(m + n),
      xout = p[l]
    )$y
  }
  q
}

#' @rdname gi
#' @export
#' @importFrom stats runif
rgi <- function(n, latent, infectious) {
  stop_if_not(
    is.numeric(n),
    length(n) == 1L,
    n >= 0,
    m = "`n` must be a non-negative number."
  )
  n <- floor(n)
  if (n == 0) {
    return(numeric(0L))
  }
  check_probs(latent)
  check_probs(infectious)
  latent <- latent / sum(latent)
  infectious <- infectious / sum(infectious)

  ## Latent period is integer-valued,
  ## equal to `i` with probability `latent[i]`
  rlat <- sample(seq_along(latent), size = n, replace = TRUE, prob = latent)

  ## Infectious waiting time is real-valued,
  ## with distribution supported on [0,length(infectious))
  ## and density constant on subintervals [i,i+1)
  rwait_floor <- sample(seq_along(infectious) - 1L, size = n, replace = TRUE,
                        prob = diwt(seq_along(infectious) - 1L, infectious = infectious))
  rwait <- runif(n, min = rwait_floor, max = rwait_floor + 1L)

  ## Sum of latent period and infectious waiting time
  ## yields generation interval
  rlat + rwait
}

## Infectious period distribution function
pip <- function(q, infectious) {
  n <- length(infectious)
  p <- q
  p[] <- NA
  p[q < 1] <- 0
  p[q >= n] <- 1
  l <- !is.na(q) & q >= 1 & q < n
  if (any(l)) {
    p[l] <- cumsum(infectious)[floor(q[l])]
  }
  p
}

## Infectious waiting time density function
diwt <- function(x, infectious) {
  d <- x
  d[] <- NA
  d[x < 0] <- 0
  l <- !is.na(x) & x >= 0
  if (any(l)) {
    d[l] <- (1 - pip(x[l], infectious = infectious)) / sum(seq_along(infectious) * infectious)
  }
  d
}
