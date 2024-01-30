library(epigrowthfit)
library(tools)
options(warn = 2L, error = if (interactive()) recover)


## finalsize ###########################################################

n <- 10L
R0 <- rlnorm(n, 0, 2)
S0 <- runif(n, 0, 1)
I0 <- runif(n, 0, 1 - S0)
Z <- S0 + emdbook::lambertW(-R0 * S0 * exp(-R0 * (S0 + I0))) / R0

stopifnot(exprs = {
	all.equal(finalsize(R0 = R0, S0 = S0, I0 = I0), Z)
	all.equal(finalsize(R0 = 0, S0 = S0, I0 = I0), double(n))
	all.equal(finalsize(R0 = Inf, S0 = S0, I0 = I0), S0)
})
assertWarning(finalsize(R0 = -1, S0 = 1, I0 = 0))
assertWarning(finalsize(R0 = 1, S0 = 1, I0 = 0.1))


## R0 ##################################################################

r <- rlnorm(10L, -3, 1)
breaks <- 0:20
probs <- diff(pgamma(breaks, shape = 1, scale = 2.5))
probs <- probs / sum(probs)

n <- length(breaks)
f <- function(r) r / sum(probs * (exp(-r * breaks[-n]) - exp(-r * breaks[-1L])) / (breaks[-1L] - breaks[-n]))

R0.1 <- function(x) R0(r = x, breaks = breaks, probs = probs)
R0.3 <- function(x) R0(r = r, breaks = breaks, probs = x)

stopifnot(exprs = {
	all.equal(R0.1(r), vapply(r, f, 0))
	all.equal(R0.1(c(0, NaN, Inf)), c(1, NaN, Inf))
	all.equal(R0.3(probs), R0.3(100 * probs))
})
assertWarning(R0.1(-1))


## timescale ###########################################################

r <- c(rlnorm(10L, -3, 1), 0, NaN, Inf)
per <- 1L
tdoubling <- timescale(r = r, per = per)
vv <- withVisible(print(tdoubling))

stopifnot(exprs = {
	identical(tdoubling, structure(log(2) / r, per = per, class = "tdoubling"))
	identical(vv[["value"]], tdoubling)
	identical(vv[["visible"]], FALSE)
})
assertWarning(timescale(-1))
