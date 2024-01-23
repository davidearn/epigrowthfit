library(epigrowthfit)
options(warn = 2L, error = if (interactive()) recover, egf.cores = 2L)


## exponential #########################################################

r <- log(2) / 20
c0 <- 100
s <- 0.2

mu <- log(c(r, c0))
Sigma <- diag(rep.int(s^2, length(mu)))

zz <- simulate(egf_model(curve = "exponential", family = "pois"),
               nsim = 20L,
               seed = 775494L,
               mu = mu,
               Sigma = Sigma,
               cstart = 10)
mm <- egf(zz,
          formula_priors = list(Sigma ~ LKJ(eta = 2)))

p1 <- as.list(coef(zz))
p2 <- as.list(coef(mm))

stopifnot(exprs = {
	max(abs(mm[["gradient"]])) < 5e-05
	all.equal(p1[["beta"]], p2[["beta"]], tolerance = 5e-02)
	all.equal(theta2cov(p1[["theta"]]), theta2cov(p2[["theta"]]), tolerance = 5e-02)
})


## subexponential ######################################################

alpha <- log(2) / 20
c0 <- 100
p <- 0.95
s <- 0.2

mu <- c(log(alpha), log(c0), qlogis(p))
Sigma <- diag(rep.int(s^2, length(mu)))

zz <- simulate(egf_model(curve = "subexponential", family = "pois"),
               nsim = 20L,
               seed = 653927L,
               mu = mu,
               Sigma = Sigma,
               cstart = 10)
mm <- egf(zz,
          formula_priors = list(beta[3L] ~ Normal(mu = qlogis(p), sigma = 0.05),
                                theta[3L] ~ Normal(mu = log(s), sigma = 0.25),
                                Sigma ~ LKJ(eta = 2)))

p1 <- as.list(coef(zz))
p2 <- as.list(coef(mm))

stopifnot(exprs = {
	max(abs(mm[["gradient"]])) < 5e-04
	all.equal(p1[["beta"]], p2[["beta"]], tolerance = 5e-02)
	all.equal(theta2cov(p1[["theta"]]), theta2cov(p2[["theta"]]), tolerance = 2e-02)
})


## gompertz ############################################################

alpha <- log(2) / 20
tinfl <- 100
K <- 25000
s <- 0.2

mu <- log(c(alpha, tinfl, K))
Sigma <- diag(rep.int(s^2, length(mu)))

zz <- simulate(egf_model(curve = "gompertz", family = "pois"),
               nsim = 20L,
               seed = 685399L,
               mu = mu,
               Sigma = Sigma,
               cstart = 10)
oo <- options(warn = 1L) # FIXME: diagnose NA/NaN function evaluation
mm <- egf(zz,
          formula_priors = list(Sigma ~ LKJ(eta = 2)))
options(oo)

p1 <- as.list(coef(zz))
p2 <- as.list(coef(mm))

stopifnot(exprs = {
	max(abs(mm[["gradient"]])) < 5e-04
	all.equal(p1[["beta"]], p2[["beta"]], tolerance = 5e-02)
	all.equal(theta2cov(p1[["theta"]]), theta2cov(p2[["theta"]]), tolerance = 2e-02)
})


## logistic ############################################################

r <- log(2) / 20
tinfl <- 100
K <- 25000
s <- 0.2

mu <- log(c(r, tinfl, K))
Sigma <- diag(rep.int(s^2, length(mu)))

zz <- simulate(egf_model(curve = "logistic", family = "pois"),
               nsim = 20L,
               seed = 397981L,
               mu = mu,
               Sigma = Sigma,
               cstart = 10)
mm <- egf(zz,
          formula_priors = list(Sigma ~ LKJ(eta = 2)))

p1 <- as.list(coef(zz))
p2 <- as.list(coef(mm))

stopifnot(exprs = {
	max(abs(mm[["gradient"]])) < 1e-02
	all.equal(p1[["beta"]], p2[["beta"]], tolerance = 1e-02)
	all.equal(theta2cov(p1[["theta"]]), theta2cov(p2[["theta"]]), tolerance = 2e-02)
})


## richards ############################################################

r <- log(2) / 20
tinfl <- 100
K <- 25000
a <- 1.005
s <- 0.2

mu <- log(c(r, tinfl, K, a))
Sigma <- diag(rep.int(s^2, length(mu)))

zz <- simulate(egf_model(curve = "richards", family = "pois"),
               nsim = 20L,
               seed = 949642L,
               mu = mu,
               Sigma = Sigma,
               cstart = 10)
mm <- egf(zz,
          formula_priors = list(beta[4L] ~ Normal(mu = log(a), sigma = 0.005),
                                theta[4L] ~ Normal(mu = log(s), sigma = 0.25),
                                Sigma ~ LKJ(eta = 2)))

p1 <- as.list(coef(zz))
p2 <- as.list(coef(mm))

stopifnot(exprs = {
	max(abs(mm[["gradient"]])) < 2e-02
	all.equal(p1[["beta"]], p2[["beta"]], tolerance = 5e-03)
	all.equal(theta2cov(p1[["theta"]]), theta2cov(p2[["theta"]]), tolerance = 2e-02)
})
