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
