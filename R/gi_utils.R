pinf <- function(q, inf) {
  n <- length(inf)
  p <- rep.int(NA_real_, length(q))
  p[q < 1] <- 0
  p[q >= n] <- 1
  is_q_in_1_n <- (!is.na(q) & q >= 1 & q < n)
  if (any(is_q_in_1_n)) {
    p[is_q_in_1_n] <- cumsum(inf)[floor(q[is_q_in_1_n])]
  }
  dim(p) <- dim(q)
  p
}

dwait <- function(x, inf) {
  n <- length(x)
  d <- rep.int(NA_real_, length(x))
  d[x < 0] <- 0
  is_x_geq_0 <- (!is.na(x) & x >= 0)
  if (any(is_x_geq_0)) {
    d[is_x_geq_0] <- (1 - pinf(x[is_x_geq_0], inf = inf)) / sum(seq_len(n) * inf)
  }
  dim(d) <- dim(x)
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
    m = paste0(
      "`lat` and `inf` must be finite numeric vectors\n",
      "of nonzero length."
    ),
    n = 2L
  )
  stop_if_not(
    lat >= 0,
    inf >= 0,
    m = "Elements of `lat` and `inf` must be non-negative.",
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
