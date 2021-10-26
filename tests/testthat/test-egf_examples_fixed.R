test_that("exponential", {
  r <- log(2) / 20
  c0 <- 100

  zz <- simulate(egf_model(curve = "exponential", family = "pois"),
    nsim = 1L,
    seed = 412575L,
    mu = log(c(r, c0)),
    cstart = 10
  )

  mm <- egf(zz)
  expect_equal(coef(mm, full = TRUE), coef(zz), tolerance = 5e-2)
})

test_that("subexponential", {
  alpha <- log(2) / 20
  c0 <- 100
  p <- 0.95

  zz <- simulate(egf_model(curve = "subexponential", family = "pois"),
    nsim = 1L,
    seed = 696182L,
    mu = c(log(alpha), log(c0), qlogis(p)),
    cstart = 10
  )

  mm <- egf(zz, formula_priors = list(logit(p) ~ Normal(mu = qlogis(p), sigma = 0.5)))
  expect_equal(coef(mm, full = TRUE), coef(zz), tolerance = 5e-2)
})

test_that("gompertz", {
  alpha <- log(2) / 20
  tinfl <- 100
  K <- 25000

  zz <- simulate(egf_model(curve = "gompertz", family = "pois"),
    nsim = 1L,
    seed = 720748L,
    mu = log(c(alpha, tinfl, K)),
    cstart = 10
  )

  mm <- egf(zz)
  expect_equal(coef(mm, full = TRUE), coef(zz), tolerance = 5e-2)
})

test_that("logistic", {
  r <- log(2) / 20
  tinfl <- 100
  K <- 25000

  zz <- simulate(egf_model(curve = "logistic", family = "pois"),
    nsim = 1L,
    seed = 366465L,
    mu = log(c(r, tinfl, K)),
    cstart = 10
  )

  mm <- egf(zz)
  expect_equal(coef(mm, full = TRUE), coef(zz), tolerance = 5e-2)
})

test_that("richards", {
  r <- log(2) / 20
  tinfl <- 100
  K <- 25000
  a <- 1.005

  zz <- simulate(egf_model(curve = "richards", family = "pois"),
    nsim = 1L,
    seed = 51520L,
    mu = log(c(r, tinfl, K, a)),
    cstart = 10
  )

  mm <- egf(zz, formula_priors = list(log(a) ~ Normal(mu = log(a), sigma = 0.05)))
  expect_equal(coef(mm, full = TRUE), coef(zz), tolerance = 5e-2)
})
