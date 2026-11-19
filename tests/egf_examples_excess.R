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


## invariance w.r.t. time offset #######################################

## time offset changed starting point of optimization when that was not
## set in the call (GH#10)
off <- 20777 # == julian(as.Date("2026-11-20"))
dts <- transform(zz$data_ts, time = time + off)
dwi <- transform(zz$data_windows, start = start + off, end = end + off)
mmOff <- update(mm, data_ts = dts, data_windows = dwi)
stopifnot(all.equal(coef(mmOff), coef(mm), tolerance = 1e-14))
##        identical, typically
