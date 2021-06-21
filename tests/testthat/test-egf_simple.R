library("epigrowthfit")
library("testthat")

test_that("exponential", {
  set.seed(412575L)

  r <- log(2) / 20
  c0 <- 100
  nbdisp <- 10

  time <- 0:100
  mu <- diff(c0 * exp(r * time))
  x <- c(NA, rnbinom(length(mu), mu = mu, size = nbdisp))
  mm <- egf(
    model = egf_model(curve = "exponential", family = "nbinom"),
    formula = x ~ time,
    formula_par = ~1,
    endpoints = data.frame(start = time[which.max(x > 0)], end = max(time))
  )
  expect_equal(unname(mm$best), log(c(r, c0, nbdisp)), tolerance = 5e-2)
})

test_that("logistic", {
  set.seed(366465L)

  r <- log(2) / 20
  c0 <- 100
  K <- 25000
  tinfl <- log(K / c0 - 1) / r
  nbdisp <- 10

  time <- 0:(ceiling(tinfl) + 1L)
  mu <- diff(K / (1 + exp(-r * (time - tinfl))))
  x <- c(NA, rnbinom(length(mu), mu = mu, size = nbdisp))
  mm <- egf(
    model = egf_model(curve = "logistic", family = "nbinom"),
    formula = x ~ time,
    formula_par = ~1,
    endpoints = data.frame(start = time[which.max(x > 0)], end = max(time))
  )
  expect_equal(unname(mm$best), log(c(r, tinfl, K, nbdisp)), tolerance = 2e-2)
})

test_that("richards", {
  set.seed(51520L)

  r <- log(2) / 20
  c0 <- 100
  K <- 25000
  a <- 1.005
  tinfl <- (log((K / c0)^a - 1) - log(a)) / (r * a)
  nbdisp <- 10

  time <- 0:(ceiling(tinfl) + 1L)
  mu <- diff(K / (1 + a * exp(-r * a * (time - tinfl))^(1 / a)))
  x <- c(NA, rnbinom(length(mu), mu = mu, size = nbdisp))
  mm <- egf(
    model = egf_model(curve = "richards", family = "nbinom"),
    formula = x ~ time,
    formula_par = ~1,
    endpoints = data.frame(start = time[which.max(x > 0)], end = max(time))
  )
  expect_equal(unname(mm$best), log(c(r, tinfl, K, a, nbdisp)), tolerance = 5e-2)
})

test_that("subexponential", {
  set.seed(696180L)

  alpha <- log(2) / 20
  c0 <- 100
  p <- 0.95
  nbdisp <- 10

  time <- 0:100
  mu <- diff((c0^(1 - p) + (1 - p) * alpha * time)^(1 / (1 - p)))
  x <- c(NA, rnbinom(length(mu), mu = mu, size = nbdisp))
  mm <- egf(
    model = egf_model(curve = "subexponential", family = "nbinom"),
    formula = x ~ time,
    formula_par = ~1,
    endpoints = data.frame(start = time[which.max(x > 0)], end = max(time))
  )
  expect_equal(unname(mm$best), c(log(alpha), log(c0), qlogis(p), log(nbdisp)), tolerance = 1e-1)
})

test_that("gompertz", {
  set.seed(720748L)

  alpha <- log(2) / 20
  c0 <- 100
  K <- 25000
  nbdisp <- 10

  time <- 0:100
  mu <- diff(K * (c0 / K)^(exp(-alpha * time)))
  x <- c(NA, rnbinom(length(mu), mu = mu, size = nbdisp))
  mm <- egf(
    model = egf_model(curve = "gompertz", family = "nbinom"),
    formula = x ~ time,
    formula_par = ~1,
    endpoints = data.frame(start = time[which.max(x > 0)], end = max(time))
  )
  expect_equal(unname(mm$best), log(c(alpha, c0, K, nbdisp)), tolerance = 5e-2)
})
