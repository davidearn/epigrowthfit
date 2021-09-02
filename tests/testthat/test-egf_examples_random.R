oo <- options(egf.cores = 4L)

test_that("exponential", {
  r <- log(2) / 20
  c0 <- 100
  disp <- 50

  zz <- simulate(egf_model(curve = "exponential", family = "nbinom"),
    nsim = 100L,
    seed = 775494L,
    mu = log(c(r, c0, disp)),
    Sigma = diag(rep_len(0.2^2, 3L)),
    cstart = 10
  )
  mm <- egf(zz,
    formula_priors_bottom = list(
      Sigma ~ LKJ(eta = 2)
    )
  )
  expect_lt(with(mm, max(abs(tmb_out$gr(best[nonrandom])))), 1e-2)
  pp <- split(data.frame(actual = zz$actual, fitted = mm$best),
              sub("\\[[0-9]+\\]$", "", names(mm$best)))
  expect_equal(pp$beta$fitted, pp$beta$actual, tolerance = 2e-2)
  expect_equal(pp$theta$fitted[1:3], pp$theta$actual[1:3], tolerance = 2e-2)
  expect_equal(pp$theta$fitted[-(1:3)], pp$theta$actual[-(1:3)], tolerance = 5e-1)
})

test_that("subexponential", {
  alpha <- log(2) / 20
  c0 <- 100
  p <- 0.95
  disp <- 50

  zz <- simulate(egf_model(curve = "subexponential", family = "nbinom"),
    nsim = 100L,
    seed = 653927L,
    mu = c(log(alpha), log(c0), qlogis(p), log(disp)),
    Sigma = diag(rep_len(0.2^2, 4L)),
    cstart = 10
  )
  mm <- egf(zz,
    formula_priors_bottom = list(
      beta[3L] ~ Normal(mu = qlogis(p), sigma = 0.02),
      theta[3L] ~ Normal(mu = log(0.2), sigma = 0.125),
      Sigma ~ LKJ(eta = 2)
    )
  )
  expect_lt(with(mm, max(abs(tmb_out$gr(best[nonrandom])))), 1e-2)
  pp <- split(data.frame(actual = zz$actual, fitted = mm$best),
              sub("\\[[0-9]+\\]$", "", names(zz$actual)))
  expect_equal(pp$beta$fitted, pp$beta$actual, tolerance = 2e-2)
  expect_equal(pp$theta$fitted[1:4], pp$theta$actual[1:4], tolerance = 1e-1)
  expect_equal(pp$theta$fitted[-(1:4)], pp$theta$actual[-(1:4)], tolerance = 1e-1)
})

test_that("gompertz", {
  alpha <- log(2) / 20
  tinfl <- 100
  K <- 25000
  disp <- 50

  zz <- simulate(egf_model(curve = "gompertz", family = "nbinom"),
    nsim = 100L,
    seed = 685398L,
    mu = log(c(alpha, tinfl, K, disp)),
    Sigma = diag(rep_len(0.2^2, 4L)),
    cstart = 10
  )
  mm <- egf(zz,
    formula_priors_bottom = list(
      Sigma ~ LKJ(eta = 2)
    )
  )
  expect_lt(with(mm, max(abs(tmb_out$gr(best[nonrandom])))), 1e-1)
  pp <- split(data.frame(actual = zz$actual, fitted = mm$best),
              sub("\\[[0-9]+\\]$", "", names(zz$actual)))
  expect_equal(pp$beta$fitted, pp$beta$actual, tolerance = 5e-2)
  expect_equal(pp$theta$fitted[1:4], pp$theta$actual[1:4], tolerance = 2e-2)
  expect_equal(pp$theta$fitted[-(1:4)], pp$theta$actual[-(1:4)], tolerance = 1)
})

test_that("logistic", {
  r <- log(2) / 20
  tinfl <- 100
  K <- 25000
  disp <- 50

  zz <- simulate(egf_model(curve = "logistic", family = "nbinom"),
    nsim = 100L,
    seed = 397981L,
    mu = log(c(r, tinfl, K, disp)),
    Sigma = diag(rep_len(0.2^2, 4L)),
    cstart = 10
  )
  mm <- egf(zz,
    formula_priors_bottom = list(
      Sigma ~ LKJ(eta = 2)
    )
  )
  expect_lt(with(mm, max(abs(tmb_out$gr(best[nonrandom])))), 1e-2)
  pp <- split(data.frame(actual = zz$actual, fitted = mm$best),
              sub("\\[[0-9]+\\]$", "", names(zz$actual)))
  expect_equal(pp$beta$fitted, pp$beta$actual, tolerance = 5e-3)
  expect_equal(pp$theta$fitted[1:4], pp$theta$actual[1:4], tolerance = 1e-1)
  expect_equal(pp$theta$fitted[-(1:4)], pp$theta$actual[-(1:4)], tolerance = 2e-1)
})

test_that("richards", {
  r <- log(2) / 20
  tinfl <- 100
  K <- 25000
  a <- 1.005
  disp <- 10

  zz <- simulate(egf_model(curve = "richards", family = "nbinom"),
    nsim = 100L,
    seed = 949642L,
    mu = log(c(r, tinfl, K, a, disp)),
    Sigma = diag(rep_len(0.2^2, 5L)),
    cstart = 10
  )
  mm <- egf(zz,
    formula_priors_bottom = list(
      beta[4L] ~ Normal(mu = log(a), sigma = 0.02),
      theta[4L] ~ Normal(mu = log(0.2), sigma = 0.005),
      Sigma ~ LKJ(eta = 2)
    )
  )
  expect_lt(with(mm, max(abs(tmb_out$gr(best[nonrandom])))), 1e-2)
  pp <- split(data.frame(actual = zz$actual, fitted = mm$best),
              sub("\\[[0-9]+\\]$", "", names(zz$actual)))
  expect_equal(pp$beta$fitted, pp$beta$actual, tolerance = 5e-3)
  expect_equal(pp$theta$fitted[1:5], pp$theta$actual[1:5], tolerance = 1e-1)
  expect_equal(pp$theta$fitted[-(1:5)], pp$theta$actual[-(1:5)], tolerance = 5e-1)
})

options(oo)
