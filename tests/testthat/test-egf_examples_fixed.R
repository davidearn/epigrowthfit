test_that("exponential", {
  r <- log(2) / 20
  c0 <- 100
  disp <- 50

  zz <- simulate(egf_model(curve = "exponential", family = "nbinom"),
    nsim = 1L,
    seed = 412575L,
    mu = log(c(r, c0, disp)),
    cstart = 10
  )
  mm <- egf(zz)
  expect_equal(mm$best, zz$actual, tolerance = 1e-1)
})

test_that("subexponential", {
  alpha <- log(2) / 20
  c0 <- 100
  p <- 0.95
  disp <- 50

  zz <- simulate(egf_model(curve = "subexponential", family = "nbinom"),
    nsim = 1L,
    seed = 696182L,
    mu = c(log(alpha), log(c0), qlogis(p), log(disp)),
    cstart = 10
  )
  mm <- egf(zz,
    formula_priors_bottom = list(
      beta[3L] ~ Normal(mu = qlogis(p), sigma = 0.5)
    )
  )
  expect_equal(mm$best, zz$actual, tolerance = 5e-2)
})

test_that("gompertz", {
  alpha <- log(2) / 20
  tinfl <- 100
  K <- 25000
  disp <- 50

  zz <- simulate(egf_model(curve = "gompertz", family = "nbinom"),
    nsim = 1L,
    seed = 720748L,
    mu = log(c(alpha, tinfl, K, disp)),
    cstart = 10
  )
  mm <- egf(zz)
  expect_equal(mm$best, zz$actual, tolerance = 5e-2)
})

test_that("logistic", {
  r <- log(2) / 20
  tinfl <- 100
  K <- 25000
  disp <- 50

  zz <- simulate(egf_model(curve = "logistic", family = "nbinom"),
    nsim = 1L,
    seed = 366465L,
    mu = log(c(r, tinfl, K, disp)),
    cstart = 10
  )
  mm <- egf(zz)
  expect_equal(mm$best, zz$actual, tolerance = 5e-2)
})

test_that("richards", {
  r <- log(2) / 20
  tinfl <- 100
  K <- 25000
  a <- 1.005
  disp <- 50

  zz <- simulate(egf_model(curve = "richards", family = "nbinom"),
    nsim = 1L,
    seed = 51520L,
    mu = log(c(r, tinfl, K, a, disp)),
    cstart = 10
  )
  mm <- egf(zz,
    formula_priors_bottom = list(
      beta[4L] ~ Normal(mu = log(a), sigma = 0.05)
    )
  )
  expect_equal(mm$best, zz$actual, tolerance = 5e-2)
})
