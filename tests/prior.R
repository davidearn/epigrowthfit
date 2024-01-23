library(epigrowthfit)
library(tools)
options(warn = 2L, error = if (interactive()) recover)


## Normal ##############################################################

mu <- rnorm(10L)
sigma <- rlnorm(5L)

stopifnot(identical(Normal(mu = mu, sigma = sigma),
                    structure(list(family = "norm",
                                   parameters = list(mu = mu, sigma = sigma)),
                              class = "egf_prior")))

assertError(Normal(mu = "1",        sigma = sigma))
assertError(Normal(mu = double(0L), sigma = sigma))
assertError(Normal(mu = NaN,        sigma = sigma))
assertError(Normal(mu = Inf,        sigma = sigma))
assertError(Normal(mu = mu, sigma = "1"))
assertError(Normal(mu = mu, sigma = double(0L)))
assertError(Normal(mu = mu, sigma = NaN))
assertError(Normal(mu = mu, sigma = Inf))
assertError(Normal(mu = mu, sigma =  0))
assertError(Normal(mu = mu, sigma = -1))


## LKJ #################################################################

eta <- c(0.1, 1, 10)

stopifnot(identical(LKJ(eta = eta),
                    structure(list(family = "lkj",
                                   parameters = list(eta = eta)),
                              class = "egf_prior")))

assertError(LKJ(eta = "1"))
assertError(LKJ(eta = double(0L)))
assertError(LKJ(eta = NaN))
assertError(LKJ(eta = Inf))
assertError(LKJ(eta = 0))
assertError(LKJ(eta = -1))


## (Inverse)?Wishart ###################################################

set.seed(230719L)
n <- 4L
A <- matrix(rnorm(n * n), n, n)
S <- crossprod(A)
R <- chol(S)
R1 <- R * rep(1 / diag(R, names = FALSE), each = n)

df <- 8
theta <- c(0.5 * log(diag(S, names = FALSE)), R1[upper.tri(R1)])
prior <-
	structure(list(family = NA_character_,
	               parameters = list(df = df, scale = list(theta))),
	          class = "egf_prior")

stopifnot(exprs = {
	identical(       Wishart(df = df, scale = list(S)),
	          replace(prior, "family", list(   "wishart")))
	identical(InverseWishart(df = df, scale = list(S)),
	          replace(prior, "family", list("invwishart")))
})

assertError(Wishart(df =  3, scale = S))
assertError(Wishart(df = df, scale = A))
assertError(Wishart(df = df, scale = diag(0, n)))
