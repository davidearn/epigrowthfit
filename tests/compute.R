library(epigrowthfit)
library(tools)
options(warn = 2L, error = if (interactive()) recover)


## compute_final_size ##################################################

n <- 10L
R0 <- rlnorm(n, 0, 2)
S0 <- runif(n, 0, 1)
I0 <- runif(n, 0, 1 - S0)
Z <- S0 + emdbook::lambertW(-R0 * S0 * exp(-R0 * (S0 + I0))) / R0

all.equal(compute_final_size(R0 = R0, S0 = S0, I0 = I0), Z)
all.equal(compute_final_size(R0 = 0, S0 = S0, I0 = I0), double(n))
all.equal(compute_final_size(R0 = Inf, S0 = S0, I0 = I0), S0)
assertWarning(compute_final_size(R0 = -1, S0 = 1, I0 = 0), "NA")
assertWarning(compute_final_size(R0 = 1, S0 = 1, I0 = 0.1), "NA")


## compute_R0 ##########################################################

r <- rlnorm(10L, -3, 1)
breaks <- 0:20
probs <- diff(pgamma(breaks, shape = 1, scale = 2.5))
probs <- probs / sum(probs)

n <- length(breaks)
f <- function(r) r / sum(probs * (exp(-r * breaks[-n]) - exp(-r * breaks[-1L])) / (breaks[-1L] - breaks[-n]))

.compute_R0 <- function(x) compute_R0(r = x, breaks = breaks, probs = probs)
all.equal(.compute_R0(r), vapply(r, f, 0))
all.equal(.compute_R0(c(0, NA, NaN, Inf)), c(1, NA, NaN, Inf))
assertWarning(.compute_R0(-1), "NA")
.compute_R0 <- function(x) compute_R0(r = r, breaks = breaks, probs = x)
all.equal(.compute_R0(probs), .compute_R0(100 * probs))


## compute_tdoubling ###################################################

r <- c(rlnorm(10L, -3, 1), 0, NA, NaN, Inf)
per <- 1L
tdoubling <- compute_tdoubling(r = r, per = per)
identical(tdoubling,
                 structure(log(2) / r, per = per, class = "tdoubling"))
assertWarning(compute_tdoubling(-1), "NA")

capture.output({
    expect_condition(print(tdoubling), regexp = NA)
    identical(print(tdoubling), tdoubling)
    expect_invisible(print(tdoubling))
})
