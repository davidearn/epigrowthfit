library("epigrowthfit")
library("testthat")
options(egf.cores = 4L)

test_that("exponential", {
  r <- log(2) / 20
  c0 <- 100
  nbdisp <- 50

  set.seed(775494L)
  zz <- egf_simulate(
    N = 100L,
    model = egf_model(curve = "exponential", family = "nbinom"),
    mu = log(c(r, c0, nbdisp)),
    Sigma = diag(rep_len(0.2^2, 3L)),
    cstart = 10
  )
  mm <- egf(zz)
  index <- split(seq_along(zz$actual), sub("\\[[0-9]+\\]$", "", names(zz$actual)))
  expect_equal(mm$best[index$beta], zz$actual[index$beta], tolerance = 2e-2)
  expect_equal(mm$best[index$theta[1:3]], zz$actual[index$theta[1:3]], tolerance = 1e-1)
  expect_equal(mm$best[index$theta[-(1:3)]], zz$actual[index$theta[-(1:3)]], tolerance = 5e-1)
})

test_that("subexponential", {
  alpha <- log(2) / 20
  c0 <- 100
  p <- 0.95
  nbdisp <- 50

  set.seed(653927L)
  zz <- egf_simulate(
    N = 100L,
    model = egf_model(curve = "subexponential", family = "nbinom"),
    mu = c(log(alpha), log(c0), qlogis(p), log(nbdisp)),
    Sigma = diag(rep_len(0.2^2, 4L)),
    cstart = 10
  )
  mm <- egf(zz)
  index <- split(seq_along(zz$actual), sub("\\[[0-9]+\\]$", "", names(zz$actual)))
  expect_equal(mm$best[index$beta], zz$actual[index$beta], tolerance = 5e-2)
  expect_equal(mm$best[index$theta[1:4]], zz$actual[index$theta[1:4]], tolerance = 5e-1)
  expect_equal(mm$best[index$theta[-(1:4)]], zz$actual[index$theta[-(1:4)]], tolerance = 1)
})

test_that("gompertz", {
  alpha <- log(2) / 20
  c0 <- 100
  K <- 25000
  nbdisp <- 50

  set.seed(685398L)
  zz <- egf_simulate(
    N = 100L,
    model = egf_model(curve = "gompertz", family = "nbinom"),
    mu = log(c(alpha, c0, K, nbdisp)),
    Sigma = diag(rep_len(0.2^2, 4L)),
    cstart = 10
  )
  mm <- egf(zz)
  index <- split(seq_along(zz$actual), sub("\\[[0-9]+\\]$", "", names(mm$best)))
  expect_equal(mm$best[index$beta], zz$actual[index$beta], tolerance = 2e-2)
  expect_equal(mm$best[index$theta[1:4]], zz$actual[index$theta[1:4]], tolerance = 5e-2)
  expect_equal(mm$best[index$theta[-(1:4)]], zz$actual[index$theta[-(1:4)]], tolerance = 2e-1)
})


test_that("logistic", {
  r <- log(2) / 20
  tinfl <- 100
  K <- 25000
  nbdisp <- 50

  set.seed(397981L)
  zz <- egf_simulate(
    N = 100L,
    model = egf_model(curve = "logistic", family = "nbinom"),
    mu = log(c(r, tinfl, K, nbdisp)),
    Sigma = diag(rep_len(0.2^2, 4L)),
    cstart = 10
  )
  mm <- egf(zz)
  index <- split(seq_along(zz$actual), sub("\\[[0-9]+\\]$", "", names(mm$best)))
  expect_equal(mm$best[index$beta], zz$actual[index$beta], tolerance = 5e-3)
  expect_equal(mm$best[index$theta[1:4]], zz$actual[index$theta[1:4]], tolerance = 1e-1)
  expect_equal(mm$best[index$theta[-(1:4)]], zz$actual[index$theta[-(1:4)]], tolerance = 1e-1)
})

test_that("richards", {
  r <- log(2) / 20
  tinfl <- 100
  K <- 25000
  a <- 1.005
  nbdisp <- 10

  set.seed(949642L)
  zz <- egf_simulate(
    N = 100L,
    model = egf_model(curve = "richards", family = "nbinom"),
    mu = log(c(r, tinfl, K, a, nbdisp)),
    Sigma = diag(rep_len(0.2^2, 5L)),
    cstart = 10
  )
  mm <- egf(zz)
  index <- split(seq_along(zz$actual), sub("\\[[0-9]+\\]$", "", names(mm$best)))
  expect_equal(mm$best[index$beta], zz$actual[index$beta], tolerance = 5e-2)
  expect_equal(mm$best[index$theta[1:5]], zz$actual[index$theta[1:5]], tolerance = 5e-1)
  expect_equal(mm$best[index$theta[-(1:5)]], zz$actual[index$theta[-(1:5)]], tolerance = 2e+1)
})
