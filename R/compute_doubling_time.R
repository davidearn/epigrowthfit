compute_doubling_time <- function(r) {
  UseMethod("compute_doubling_time", r)
}

compute_doubling_time.default <- function(r) {
  check(r,
    what = "numeric",
    val = c(0, NA),
    "`r` must be a numeric vector with non-negative elements."
  )
  structure(log(2) / r, class = c("doubling_time", "numeric"))
}
