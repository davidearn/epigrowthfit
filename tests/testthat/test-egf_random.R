library("epigrowthfit")
library("testthat")

test_that("exponential", {
  set.seed(775494L)

  r <- log(2) / 20
  c0 <- 100
  nbdisp <- 10

  beta <- log(c(r, c0, nbdisp))
  n <- length(beta)
  sd_b <- rep_len(0.1, n)
  theta <- c(log(sd_b), rep_len(0, n * (n - 1) / 2))
  N <- 500L
  b <- rnorm(n * N)
  Y <- matrix(beta + sd_b * b, ncol = n, byrow = TRUE)
  # colnames(Y) <- c("log(r)", "log(c0)", "log(nbdisp)")

  fp <- function(time, par) {
    r <- exp(par[1L])
    c0 <- exp(par[2L])
    diff(c0 * exp(r * time))
  }
  fs <- function(mu, size) {
    c(NA, rnbinom(length(mu), mu = mu, size = size))
  }

  time <- 0:100
  mu <- apply(Y, 1L, fp, time = time, simplify = FALSE)
  x <- Map(fs, mu = mu, size = exp(Y[, 3L]))
  dd <- data.frame(
    time = time,
    x = unlist(x, FALSE, FALSE),
    ts = gl(N, length(time))
  )
  ep <- data.frame(
    start = mapply(`[[`, list(time), vapply(x, function(n) which.max(n > 0), 0L)),
    end = max(time),
    ts = gl(N, 1L)
  )
  mm <- egf(
    model = egf_model(curve = "exponential", family = "nbinom"),
    formula = x ~ time | ts,
    formula_par = ~(1 | ts),
    data = dd,
    data_par = ep,
    endpoints = ep,
    control = egf_control(
      inner_optimizer = egf_inner_optimizer(args = list(trace = FALSE)),
      omp_num_threads = 4L
    )
  )
  pp <- split(unname(mm$best), sub("\\[[0-9]+\\]$", "", names(mm$best)))

  ## Means
  expect_equal(pp$beta, beta, tolerance = 2e-2)
  ## Log standard deviations
  expect_equal(pp$theta[1:n], theta[1:n], tolerance = 1e-2)
  ## Other covariance parameters
  expect_equal(pp$theta[-(1:n)], theta[-(1:n)], tolerance = 5e-2)

  ## Predicted values
  Y_hat <- mm$tmb_out$report(mm$best)$Y
  mu_hat <- apply(Y_hat, 1L, fp, time = time, simplify = FALSE)
  expect_equal(unlist(mu_hat, FALSE, FALSE), unlist(mu, FALSE, FALSE), tolerance = 1e-1)
})

test_that("logistic", {
  set.seed(397981L)

  r <- log(2) / 20
  c0 <- 100
  K <- 25000
  tinfl <- log(K / c0 - 1) / r
  nbdisp <- 10

  beta <- log(c(r, tinfl, K, nbdisp))
  n <- length(beta)
  sd_b <- rep_len(0.1, n)
  theta <- c(log(sd_b), rep_len(0, n * (n - 1) / 2))
  N <- 500L
  b <- rnorm(n * N)
  Y <- matrix(beta + sd_b * b, ncol = n, byrow = TRUE)
  # colnames(Y) <- c("log(r)", "log(tinfl)", "log(K)", "log(nbdisp)")

  fp <- function(time, par) {
    r <- exp(par[1L])
    tinfl <- exp(par[2L])
    K <- exp(par[3L])
    diff(K / (1 + exp(-r * (time - tinfl))))
  }
  fs <- function(mu, size) {
    c(NA, rnbinom(length(mu), mu = mu, size = size))
  }

  time <- Map(seq.int, from = 0, to = ceiling(exp(Y[, 2L])) + 1)
  mu <- Map(fp, time = time, par = as.data.frame(t(Y)))
  x <- Map(fs, mu = mu, size = exp(Y[, 4L]))
  dd <- data.frame(
    time = unlist(time, FALSE, FALSE),
    x = unlist(x, FALSE, FALSE),
    ts = rep.int(gl(N, 1L), lengths(time))
  )
  ep <- data.frame(
    start = mapply(`[[`, time, vapply(x, function(n) which.max(n > 0), 0L)),
    end = vapply(time, max, 0),
    ts = gl(N, 1L)
  )
  mm <- egf(
    model = egf_model(curve = "logistic", family = "nbinom"),
    formula = x ~ time | ts,
    formula_par = ~(1 | ts),
    data = dd,
    data_par = ep,
    endpoints = ep,
    control = egf_control(
      inner_optimizer = egf_inner_optimizer(args = list(trace = FALSE)),
      omp_num_threads = 4L
    )
  )
  pp <- split(unname(mm$best), sub("\\[[0-9]+\\]$", "", names(mm$best)))

  ## Means
  expect_equal(pp$beta, beta, tolerance = 1e-2)
  ## Log standard deviations
  expect_equal(pp$theta[1:n], theta[1:n], tolerance = 1e-1)
  ## Other covariance parameters
  expect_equal(pp$theta[-(1:n)], theta[-(1:n)], tolerance = 1e-1)

  ## Predicted values
  Y_hat <- mm$tmb_out$report(mm$best)$Y
  mu_hat <- Map(fp, time = time, par = as.data.frame(t(Y_hat)))
  expect_equal(unlist(mu_hat, FALSE, FALSE), unlist(mu, FALSE, FALSE), tolerance = 5e-1)
})

test_that("richards", {
  set.seed(949642L)

  r <- log(2) / 20
  c0 <- 100
  K <- 25000
  a <- 1.005
  tinfl <- (log((K / c0)^a - 1) - log(a)) / (r * a)
  nbdisp <- 10

  beta <- log(c(r, tinfl, K, a, nbdisp))
  n <- length(beta)
  sd_b <- rep_len(0.1, n)
  theta <- c(log(sd_b), rep_len(0, n * (n - 1) / 2))
  N <- 500L
  b <- rnorm(n * N)
  Y <- matrix(beta + sd_b * b, ncol = n, byrow = TRUE)
  # colnames(Y) <- c("log(r)", "log(tinfl)", "log(K)", "log(a)", "log(nbdisp)")

  fp <- function(time, par) {
    r <- exp(par[1L])
    tinfl <- exp(par[2L])
    K <- exp(par[3L])
    a <- exp(par[4L])
    diff(K / (1 + a * exp(-r * a * (time - tinfl))^(1 / a)))
  }
  fs <- function(mu, size) {
    c(NA, rnbinom(length(mu), mu = mu, size = size))
  }

  time <- Map(seq.int, from = 0, to = ceiling(exp(Y[, 2L])) + 1)
  mu <- Map(fp, time = time, par = as.data.frame(t(Y)))
  x <- Map(fs, mu = mu, size = exp(Y[, 5L]))
  dd <- data.frame(
    time = unlist(time, FALSE, FALSE),
    x = unlist(x, FALSE, FALSE),
    ts = rep.int(gl(N, 1L), lengths(time))
  )
  ep <- data.frame(
    start = mapply(`[[`, time, vapply(x, function(n) which.max(n > 0), 0L)),
    end = vapply(time, max, 0),
    ts = gl(N, 1L)
  )
  mm <- egf(
    model = egf_model(curve = "richards", family = "nbinom"),
    formula = x ~ time | ts,
    formula_par = ~(1 | ts),
    data = dd,
    data_par = ep,
    endpoints = ep,
    control = egf_control(
      inner_optimizer = egf_inner_optimizer(args = list(trace = FALSE)),
      omp_num_threads = 4L
    )
  )
  pp <- split(unname(mm$best), sub("\\[[0-9]+\\]$", "", names(mm$best)))

  ## Means
  expect_equal(pp$beta, beta, tolerance = 5e-2)
  ## Log standard deviations
  # expect_equal(pp$theta[1:n], theta[1:n], tolerance = 1)
  ## Other covariance parameters
  # expect_equal(pp$theta[-(1:n)], theta[-(1:n)], tolerance = 5e-1)

  ## Predicted values
  Y_hat <- mm$tmb_out$report(mm$best)$Y
  mu_hat <- Map(fp, time = time, par = as.data.frame(t(Y_hat)))
  expect_equal(unlist(mu_hat, FALSE, FALSE), unlist(mu, FALSE, FALSE), tolerance = 5e-1)
})

test_that("subexponential", {
  set.seed(653927L)

  alpha <- log(2) / 20
  c0 <- 100
  p <- 0.95
  nbdisp <- 10

  beta <- c(log(alpha), log(c0), qlogis(p), log(nbdisp))
  n <- length(beta)
  sd_b <- rep_len(0.1, n)
  theta <- c(log(sd_b), rep_len(0, n * (n - 1) / 2))
  N <- 500L
  b <- rnorm(n * N)
  Y <- matrix(beta + sd_b * b, ncol = n, byrow = TRUE)
  # colnames(Y) <- c("log(alpha)", "log(c0)", "logit(p)", "log(nbdisp)")

  fp <- function(time, par) {
    r <- exp(par[1L])
    c0 <- exp(par[2L])
    diff((c0^(1 - p) + (1 - p) * alpha * time)^(1 / (1 - p)))
  }
  fs <- function(mu, size) {
    c(NA, rnbinom(length(mu), mu = mu, size = size))
  }

  time <- 0:100
  mu <- apply(Y, 1L, fp, time = time, simplify = FALSE)
  x <- Map(fs, mu = mu, size = exp(Y[, 4L]))
  dd <- data.frame(
    time = time,
    x = unlist(x, FALSE, FALSE),
    ts = gl(N, length(time))
  )
  ep <- data.frame(
    start = mapply(`[[`, list(time), vapply(x, function(n) which.max(n > 0), 0L)),
    end = max(time),
    ts = gl(N, 1L)
  )
  mm <- egf(
    model = egf_model(curve = "subexponential", family = "nbinom"),
    formula = x ~ time | ts,
    formula_par = ~(1 | ts),
    data = dd,
    data_par = ep,
    endpoints = ep,
    control = egf_control(
      inner_optimizer = egf_inner_optimizer(args = list(trace = FALSE)),
      omp_num_threads = 4L
    )
  )
  pp <- split(unname(mm$best), sub("\\[[0-9]+\\]$", "", names(mm$best)))

  ## Means
  expect_equal(pp$beta, beta, tolerance = 5e-2)
  ## Log standard deviations
  expect_equal(pp$theta[1:n], theta[1:n], tolerance = 5e-1)
  ## Other covariance parameters
  # expect_equal(pp$theta[-(1:n)], theta[-(1:n)], tolerance = 5)

  ## Predicted values
  Y_hat <- mm$tmb_out$report(mm$best)$Y
  mu_hat <- apply(Y_hat, 1L, fp, time = time, simplify = FALSE)
  expect_equal(unlist(mu_hat, FALSE, FALSE), unlist(mu, FALSE, FALSE), tolerance = 1e-1)
})

test_that("gompertz", {
  set.seed(685398L)

  alpha <- log(2) / 20
  c0 <- 100
  K <- 25000
  nbdisp <- 10

  beta <- log(c(alpha, c0, K, nbdisp))
  n <- length(beta)
  sd_b <- rep_len(0.1, n)
  theta <- c(log(sd_b), rep_len(0, n * (n - 1) / 2))
  N <- 500L
  b <- rnorm(n * N)
  Y <- matrix(beta + sd_b * b, ncol = n, byrow = TRUE)
  # colnames(Y) <- c("log(alpha)", "log(c0)", "log(K)", "log(nbdisp)")

  fp <- function(time, par) {
    alpha <- exp(par[1L])
    c0 <- exp(par[2L])
    K <- exp(par[3L])
    diff(K * (c0 / K)^(exp(-alpha * time)))
  }
  fs <- function(mu, size) {
    c(NA, rnbinom(length(mu), mu = mu, size = size))
  }

  time <- Map(seq.int, from = 0, to = ceiling((1 / exp(Y[, 1L])) * log(Y[, 3L] - Y[, 2L])) + 1)
  mu <- Map(fp, time = time, par = as.data.frame(t(Y)))
  x <- Map(fs, mu = mu, size = exp(Y[, 4L]))
  dd <- data.frame(
    time = unlist(time, FALSE, FALSE),
    x = unlist(x, FALSE, FALSE),
    ts = rep.int(gl(N, 1L), lengths(time))
  )
  ep <- data.frame(
    start = mapply(`[[`, time, vapply(x, function(n) which.max(n > 0), 0L)),
    end = vapply(time, max, 0),
    ts = gl(N, 1L)
  )
  mm <- egf(
    model = egf_model(curve = "gompertz", family = "nbinom"),
    formula = x ~ time | ts,
    formula_par = ~(1 | ts),
    data = dd,
    data_par = ep,
    endpoints = ep,
    control = egf_control(
      inner_optimizer = egf_inner_optimizer(args = list(trace = FALSE)),
      omp_num_threads = 4L
    )
  )
  pp <- split(unname(mm$best), sub("\\[[0-9]+\\]$", "", names(mm$best)))

  ## Means
  expect_equal(pp$beta, beta, tolerance = 5e-2)
  ## Log standard deviations
  expect_equal(pp$theta[1:n], theta[1:n], tolerance = 1e-01)
  ## Other covariance parameters
  expect_equal(pp$theta[-(1:n)], theta[-(1:n)], tolerance = 2e-1)

  ## Predicted values
  Y_hat <- mm$tmb_out$report(mm$best)$Y
  mu_hat <- Map(fp, time = time, par = as.data.frame(t(Y_hat)))
  expect_equal(unlist(mu_hat, FALSE, FALSE), unlist(mu, FALSE, FALSE), tolerance = 1e-1)
})
