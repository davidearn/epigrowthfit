test_that("day_of_week", {
  r <- log(2) / 20
  tinfl <- 100
  K <- 25000
  w <- c(1, 1, 1, 0.8, 0.8, 1.4)
  disp <- 50

  zz <- simulate(egf_model(curve = "logistic", family = "nbinom", day_of_week = 3L),
    nsim = 1L,
    seed = 366465L,
    mu = log(c(r, tinfl, K, disp, w)),
    cstart = 10
  )
  mm <- egf(zz,
    formula_priors = list(
      beta[4:9] ~ Normal(mu = log(w), sigma = 0.2)
    )
  )
  expect_equal(mm$best, zz$actual, tolerance = 2e-1, ignore_attr = "lengths")
})
