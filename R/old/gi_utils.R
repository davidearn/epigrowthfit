pinf <- function(q, inf) {
  n <- length(inf)
  p <- q
  p[] <- NA
  p[q < 1] <- 0
  p[q >= n] <- 1
  l <- !is.na(q) & q >= 1 & q < n
  if (any(l)) {
    p[l] <- cumsum(inf)[floor(q[l])]
  }
  p
}

dwait <- function(x, inf) {
  d <- x
  d[] <- NA
  d[x < 0] <- 0
  l <- !is.na(x) & x >= 0
  if (any(l)) {
    d[l] <- (1 - pinf(x[l], inf = inf)) / sum(seq_along(inf) * inf)
  }
  d
}

check_gi <- function(lat, inf) {
  stop_if_not(
    is.numeric(lat),
    is.numeric(inf),
    length(lat) > 0L,
    length(inf) > 0L,
    is.finite(lat),
    is.finite(inf),
    m = "`lat` and `inf` must be finite numeric vectors\nof nonzero length.",
    n = 2L
  )
  stop_if_not(
    lat >= 0,
    inf >= 0,
    m = "`lat` and `inf` must be non-negative.",
    n = 2L
  )
  stop_if_not(
    any(lat > 0),
    any(inf > 0),
    m = "`lat` and `inf` must each sum to a positive number.",
    n = 2L
  )
  invisible(NULL)
}
