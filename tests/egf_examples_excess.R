library(epigrowthfit)
options(warn = 2L, error = if (interactive()) recover)


## excess ##############################################################

r <- log(2) / 20
tinfl <- 100
K <- 25000
b <- 10
disp <- 50

zz <- simulate(egf_model(curve = "logistic", family = "nbinom",
                         excess = TRUE),
               nsim = 1L,
               seed = 366465L,
               mu = log(c(r, tinfl, K, b, disp)),
               cstart = 10)
mm <- egf(zz, formula_priors = list(log(b) ~ Normal(mu = 2.5, sigma = 1)))

stopifnot(all.equal(coef(zz), coef(mm), tolerance = 5e-02))
