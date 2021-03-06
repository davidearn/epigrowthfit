library("epigrowthfit")
library("testthat")

test_that("exponential", {
  r <- log(2) / 20
  c0 <- 100
  nbdisp <- 50

  set.seed(412575L)
  zz <- egf_simulate(
    N = 1L,
    model = egf_model(curve = "exponential", family = "nbinom"),
    mu = log(c(r, c0, nbdisp)),
    cstart = 10
  )
  mm <- egf(zz)
  expect_equal(mm$best, zz$actual, tolerance = 1e-1)
})

test_that("subexponential", {
  alpha <- log(2) / 20
  c0 <- 100
  p <- 0.95
  nbdisp <- 50

  set.seed(696182L)
  zz <- egf_simulate(
    N = 1L,
    model = egf_model(curve = "subexponential", family = "nbinom"),
    mu = c(log(alpha), log(c0), qlogis(p), log(nbdisp)),
    cstart = 10
  )
  mm <- egf(zz)
  expect_equal(mm$best, zz$actual, tolerance = 5e-1)
})

test_that("gompertz", {
  alpha <- log(2) / 20
  c0 <- 100
  K <- 25000
  nbdisp <- 50

  set.seed(720748L)
  zz <- egf_simulate(
    N = 1L,
    model = egf_model(curve = "gompertz", family = "nbinom"),
    mu = log(c(alpha, c0, K, nbdisp)),
    cstart = 10
  )
  mm <- egf(zz)
  expect_equal(mm$best, zz$actual, tolerance = 1e-2)
})

test_that("logistic", {
  r <- log(2) / 20
  tinfl <- 100
  K <- 25000
  nbdisp <- 50

  set.seed(366465L)
  zz <- egf_simulate(
    N = 1L,
    model = egf_model(curve = "logistic", family = "nbinom"),
    mu = log(c(r, tinfl, K, nbdisp)),
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
  nbdisp <- 50

  set.seed(51520L)
  zz <- egf_simulate(
    N = 1L,
    model = egf_model(curve = "richards", family = "nbinom"),
    mu = log(c(r, tinfl, K, a, nbdisp)),
    cstart = 10
  )
  mm <- egf(zz)
  expect_equal(mm$best, zz$actual, tolerance = 1e-1)
})
