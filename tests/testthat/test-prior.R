test_that("Normal", {
  mu <- rnorm(10L)
  sigma <- rlnorm(5L)
  prior <- Normal(mu = mu, sigma = sigma)
  expect_type(prior, "list")
  expect_s3_class(prior, "egf_prior")
  expect_identical(unclass(prior), list(family = "norm", parameters = list(mu = mu, sigma = sigma)))

  expect_error(Normal(mu = "1",       sigma = sigma))
  expect_error(Normal(mu = numeric(0L), sigma = sigma))
  expect_error(Normal(mu = NA_real_,    sigma = sigma))
  expect_error(Normal(mu = Inf,         sigma = sigma))

  expect_error(Normal(mu = mu, sigma = "1"))
  expect_error(Normal(mu = mu, sigma = numeric(0L)))
  expect_error(Normal(mu = mu, sigma = NA_real_))
  expect_error(Normal(mu = mu, sigma = Inf))
  expect_error(Normal(mu = mu, sigma = 0))
  expect_error(Normal(mu = mu, sigma = -1))
})

test_that("LKJ", {
  eta <- c(0.1, 1, 10)
  prior <- LKJ(eta = eta)
  expect_type(prior, "list")
  expect_s3_class(prior, "egf_prior")
  expect_identical(unclass(prior), list(family = "lkj", parameters = list(eta = eta)))

  expect_error(LKJ(eta = "1"))
  expect_error(LKJ(eta = numeric(0L)))
  expect_error(LKJ(eta = NA_real_))
  expect_error(LKJ(eta = Inf))
  expect_error(LKJ(eta = 0))
  expect_error(LKJ(eta = -1))
})

test_that("(Inverse)?Wishart", {
  df <- 8
  A <- matrix(rnorm(16L), 4L, 4L)
  scale <- A %*% t(A)
  prior <- Wishart(df = df, scale = scale)

  log_sd <- 0.5 * log(diag(scale))
  R <- chol(scale)
  R[] <- R * rep(1 / diag(R), each = nrow(R))
  chol <- R[upper.tri(R)]

  expect_type(prior, "list")
  expect_s3_class(prior, "egf_prior")
  expect_identical(unclass(prior), list(family = "wishart", parameters = list(df = df, scale = list(c(log_sd, chol)))))

  expect_error(Wishart(df = 3, scale = scale)) # 'df' not greater than 'n - 1'
  expect_error(Wishart(df = df, scale = A)) # 'scale' not symmetric
  expect_error(Wishart(df = df, scale = diag(0:3))) # 'scale' not positive definite

  expect_identical(Wishart(df = df, scale = list(scale)), prior)
  expect_identical(InverseWishart(df = df, scale = scale), replace(prior, "family", list("invwishart")))
})
