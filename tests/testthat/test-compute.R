test_that("compute_final_size", {
  n <- 10L
  R0 <- rlnorm(n, 0, 2)
  S0 <- runif(n, 0, 1)
  I0 <- runif(n, 0, 1 - S0)
  Z <- S0 + emdbook::lambertW(-R0 * S0 * exp(-R0 * (S0 + I0))) / R0

  expect_equal(compute_final_size(R0 = R0, S0 = S0, I0 = I0), Z)
  expect_equal(compute_final_size(R0 = 0, S0 = S0, I0 = I0), rep_len(0, n))
  expect_equal(compute_final_size(R0 = Inf, S0 = S0, I0 = I0), S0)
  expect_warning(compute_final_size(R0 = -1, S0 = 1, I0 = 0), "NA")
  expect_warning(compute_final_size(R0 = 1, S0 = 1, I0 = 0.1), "NA")
})

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

test_that("compute_tdoubling", {
  r <- c(rlnorm(10L, -3, 1), 0, NA, NaN, Inf)
  per <- 1L
  tdoubling <- compute_tdoubling(r = r, per = per)

  expect_equal(as.numeric(tdoubling), log(2) / r)
  expect_s3_class(tdoubling, "tdoubling")
  expect_equal(attr(tdoubling, "per"), per)
  expect_warning(compute_tdoubling(-1), "NA")

  expect_equal(print(tdoubling), as.numeric(tdoubling))
  expect_invisible(print(tdoubling))
})
