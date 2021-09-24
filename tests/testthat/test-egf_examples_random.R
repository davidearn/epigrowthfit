oo <- options(egf.cores = 4L)

test_that("exponential", {
  r <- log(2) / 20
  c0 <- 100

  mu <- log(c(r, c0))
  Sigma <- diag(rep_len(0.2^2, length(mu)))

  zz <- simulate(egf_model(curve = "exponential", family = "pois"),
    nsim = 20L,
    seed = 775494L,
    mu = mu,
    Sigma = Sigma,
    cstart = 10
  )
  mm <- egf(zz,
    formula_priors_bottom = list(
      Sigma ~ LKJ(eta = 2)
    )
  )

  expect_lt(max(abs(mm$gradient)), 5e-5)
  pp <- split(data.frame(actual = zz$actual, fitted = mm$best),
              sub("\\[[0-9]+\\]$", "", names(mm$best)))
  expect_equal(pp$beta$fitted, pp$beta$actual, tolerance = 5e-2)
  expect_equal(pp$theta$fitted[seq_along(mu)], pp$theta$actual[seq_along(mu)], tolerance = 2e-1)
  expect_equal(pp$theta$fitted[-seq_along(mu)], pp$theta$actual[-seq_along(mu)], tolerance = 5e-1)
})

test_that("subexponential", {
  alpha <- log(2) / 20
  c0 <- 100
  p <- 0.95

  mu <- c(log(alpha), log(c0), qlogis(p))
  Sigma <- diag(rep_len(0.2^2, length(mu)))

  zz <- simulate(egf_model(curve = "subexponential", family = "pois"),
    nsim = 20L,
    seed = 653927L,
    mu = mu,
    Sigma = Sigma,
    cstart = 10
  )
  mm <- egf(zz,
    formula_priors_bottom = list(
      beta[3L] ~ Normal(mu = mu[3L], sigma = 0.05),
      theta[3L] ~ Normal(mu = 0.5 * log(Sigma[3L, 3L]), sigma = 0.25),
      Sigma ~ LKJ(eta = 2)
    )
  )

  expect_lt(max(abs(mm$gradient)), 5e-4)
  pp <- split(data.frame(actual = zz$actual, fitted = mm$best),
              sub("\\[[0-9]+\\]$", "", names(zz$actual)))
  expect_equal(pp$beta$fitted, pp$beta$actual, tolerance = 5e-2)
  expect_equal(pp$theta$fitted[seq_along(mu)], pp$theta$actual[seq_along(mu)], tolerance = 2e-1)
  expect_equal(pp$theta$fitted[-seq_along(mu)], pp$theta$actual[-seq_along(mu)], tolerance = 5e-1)
})

test_that("gompertz", {
  alpha <- log(2) / 20
  tinfl <- 100
  K <- 25000

  mu <- log(c(alpha, tinfl, K))
  Sigma <- diag(rep_len(0.2^2, length(mu)))

  zz <- simulate(egf_model(curve = "gompertz", family = "pois"),
    nsim = 20L,
    seed = 685399L,
    mu = mu,
    Sigma = Sigma,
    cstart = 10
  )
  mm <- egf(zz,
    formula_priors_bottom = list(
      Sigma ~ LKJ(eta = 2)
    )
  )

  expect_lt(max(abs(mm$gradient)), 5e-4)
  pp <- split(data.frame(actual = zz$actual, fitted = mm$best),
              sub("\\[[0-9]+\\]$", "", names(zz$actual)))
  expect_equal(pp$beta$fitted, pp$beta$actual, tolerance = 5e-2)
  expect_equal(pp$theta$fitted[seq_along(mu)], pp$theta$actual[seq_along(mu)], tolerance = 2e-1)
  expect_equal(pp$theta$fitted[-seq_along(mu)], pp$theta$actual[-seq_along(mu)], tolerance = 2)
})

test_that("logistic", {
  r <- log(2) / 20
  tinfl <- 100
  K <- 25000

  mu <- log(c(r, tinfl, K))
  Sigma <- diag(rep_len(0.2^2, length(mu)))

  zz <- simulate(egf_model(curve = "logistic", family = "pois"),
    nsim = 20L,
    seed = 397981L,
    mu = mu,
    Sigma = Sigma,
    cstart = 10
  )
  mm <- egf(zz,
    formula_priors_bottom = list(
      Sigma ~ LKJ(eta = 2)
    )
  )

  expect_lt(max(abs(mm$gradient)), 1e-2)
  pp <- split(data.frame(actual = zz$actual, fitted = mm$best),
              sub("\\[[0-9]+\\]$", "", names(zz$actual)))
  expect_equal(pp$beta$fitted, pp$beta$actual, tolerance = 1e-2)
  expect_equal(pp$theta$fitted[seq_along(mu)], pp$theta$actual[seq_along(mu)], tolerance = 5e-2)
  expect_equal(pp$theta$fitted[-seq_along(mu)], pp$theta$actual[-seq_along(mu)], tolerance = 1e-1)
})

test_that("richards", {
  r <- log(2) / 20
  tinfl <- 100
  K <- 25000
  a <- 1.005

  mu <- log(c(r, tinfl, K, a))
  Sigma <- diag(rep_len(0.2^2, length(mu)))

  zz <- simulate(egf_model(curve = "richards", family = "pois"),
    nsim = 20L,
    seed = 949642L,
    mu = mu,
    Sigma = Sigma,
    cstart = 10
  )
  mm <- egf(zz,
    formula_priors_bottom = list(
      beta[4L] ~ Normal(mu = log(a), sigma = 0.005),
      theta[4L] ~ Normal(mu = log(0.2), sigma = 0.25),
      Sigma ~ LKJ(eta = 2)
    )
  )
  expect_lt(max(abs(mm$gradient)), 2e-2)
  pp <- split(data.frame(actual = zz$actual, fitted = mm$best),
              sub("\\[[0-9]+\\]$", "", names(zz$actual)))
  expect_equal(pp$beta$fitted, pp$beta$actual, tolerance = 5e-3)
  expect_equal(pp$theta$fitted[seq_along(mu)], pp$theta$actual[seq_along(mu)], tolerance = 1e-1)
  expect_equal(pp$theta$fitted[-seq_along(mu)], pp$theta$actual[-seq_along(mu)], tolerance = 2e-1)
})

options(oo)
