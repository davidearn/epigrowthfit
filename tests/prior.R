library(epigrowthfit)
options(warn = 2L, error = recover)


## Normal ######
    mu <- rnorm(10L)
    sigma <- rlnorm(5L)
    prior <- Normal(mu = mu, sigma = sigma)
    expect_type(prior, "list")
    expect_s3_class(prior, "egf_prior", exact = TRUE)
    identical(unclass(prior),
                     list(family = "norm",
                          parameters = list(mu = mu, sigma = sigma)))

    assertError(Normal(mu = "1",        sigma = sigma))
    assertError(Normal(mu = double(0L), sigma = sigma))
    assertError(Normal(mu = NA_real_,   sigma = sigma))
    assertError(Normal(mu = Inf,        sigma = sigma))

    assertError(Normal(mu = mu, sigma = "1"))
    assertError(Normal(mu = mu, sigma = double(0L)))
    assertError(Normal(mu = mu, sigma = NA_real_))
    assertError(Normal(mu = mu, sigma = Inf))
    assertError(Normal(mu = mu, sigma = 0))
    assertError(Normal(mu = mu, sigma = -1))


## LKJ ######
    eta <- c(0.1, 1, 10)
    prior <- LKJ(eta = eta)
    expect_type(prior, "list")
    expect_s3_class(prior, "egf_prior", exact = TRUE)
    identical(unclass(prior),
                     list(family = "lkj", parameters = list(eta = eta)))

    assertError(LKJ(eta = "1"))
    assertError(LKJ(eta = double(0L)))
    assertError(LKJ(eta = NA_real_))
    assertError(LKJ(eta = Inf))
    assertError(LKJ(eta = 0))
    assertError(LKJ(eta = -1))


## (Inverse)?Wishart ######
    df <- 8
    A <- matrix(rnorm(16L), 4L, 4L)
    scale <- A %*% t(A)
    prior <- Wishart(df = df, scale = scale)

    log_sd <- 0.5 * log(diag(scale))
    R <- chol(scale)
    R[] <- R * rep(1 / diag(R), each = nrow(R))
    chol <- R[upper.tri(R)]

    expect_type(prior, "list")
    expect_s3_class(prior, "egf_prior", exact = TRUE)
    all.equal(unclass(prior),
                 list(family = "wishart",
                      parameters = list(df = df, scale=list(c(log_sd, chol)))))

    assertError(Wishart(df = 3, scale = scale)) # 'df' not greater than 'n - 1'
    assertError(Wishart(df = df, scale = A)) # 'scale' not symmetric
    assertError(Wishart(df = df, scale = diag(0:3))) # 'scale' not positive definite

    identical(Wishart(df = df, scale = list(scale)),
                     prior)
    identical(InverseWishart(df = df, scale = scale),
                     replace(prior, "family", list("invwishart")))

