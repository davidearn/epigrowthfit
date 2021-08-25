test_that("compute_R0", {
  r <- rlnorm(10L, -3, 1)
  breaks <- 0:20
  probs <- diff(pgamma(breaks, shape = 1, scale = 2.5))
  probs <- probs / sum(probs)

  n <- length(breaks)
  f <- function(r) {
    r / sum(probs * (exp(-r * breaks[-n]) - exp(-r * breaks[-1L])) / (breaks[-1L] - breaks[-n]))
  }

  compute_R0_ <- function(x) compute_R0(r = x, breaks = breaks, probs = probs)
  expect_equal(compute_R0_(r), vapply(r, f, 0))
  expect_equal(compute_R0_(c(0, NA, NaN, Inf)), c(1, NA, NaN, Inf))
  expect_warning(compute_R0_(-1), "NA")
  compute_R0_ <- function(x) compute_R0(r = r, breaks = breaks, probs = x)
  expect_equal(compute_R0_(probs), compute_R0_(100 * probs))
})
