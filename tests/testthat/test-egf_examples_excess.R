test_that("excess", {
  r <- log(2) / 20
  tinfl <- 100
  K <- 25000
  b <- 10
  disp <- 50

  zz <- simulate(egf_model(curve = "logistic", excess = TRUE, family = "nbinom"),
    nsim = 1L,
    seed = 366465L,
    mu = log(c(r, tinfl, K, b, disp)),
    cstart = 10
  )
  mm <- egf(zz,
    formula_priors_top = list(
      log(b) ~ Normal(mu = log(b), sigma = 0.5)
    )
  )
  expect_equal(mm$best, zz$actual, tolerance = 5e-2)
})
