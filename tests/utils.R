attach(asNamespace("epigrowthfit"))
options(warn = 2L, error = if (interactive()) recover)


## disambiguate ########################################################

x <- c("c", "b", "a", "b", "a",
       "a", "b", "c", "c", "c")
y <- c("c[1]", "b[1]", "a[1]", "b[2]", "a[2]",
       "a[3]", "b[3]", "c[2]", "c[3]", "c[4]")
stopifnot(identical(disambiguate(x), y))
x <- `names<-`(seq_along(x), x)
y <- `names<-`(seq_along(x), y)
stopifnot(identical(disambiguate(x, nms = TRUE), y))


## literal_rle #########################################################

x <- c(0, NA_real_, NaN, Inf, 1)
times <- 1:5
y <- rep.int(x, times)
rle_y <- literal_rle(y)
stopifnot(exprs = {
	identical(rle_y, list(lengths = times, values = x))
	identical(y, inverse.rle(rle_y))
})


## locf ################################################################

x <- c(NA, NA, 1, NA, 2, 2, 3, NA)
stopifnot(exprs = {
	identical(locf(x), c(NA, NA, 1, 1, 2, 2, 3, 3))
	identical(locf(x, x0 = 0), c(0, 0, 1, 1, 2, 2, 3, 3))
})


## wald ################################################################

estimate <- rnorm(6L, 0, 1)
se <- rlnorm(6L, 0, 0.1)
level <- 0.95
q <- qchisq(level, df = 1)
W <- wald(estimate = estimate, se = se, level = level)
stopifnot(exprs = {
	is.double(W)
	identical(dim(W), c(6L, 2L))
	identical(dimnames(W), list(NULL, c("lower", "upper")))
	all.equal(W[, "lower"], estimate - sqrt(q) * se)
	all.equal(W[, "upper"], estimate + sqrt(q) * se)
})


## cov2theta, theta2cov ################################################

set.seed(230719L)
n <- 6L
A <- matrix(rnorm(n * n), n, n)
S <- crossprod(A)
R <- chol(S)
R1 <- R * rep(1 / diag(R, names = FALSE), each = n)

theta <- c(0.5 * log(diag(S, names = FALSE)), R1[upper.tri(R1)])
stopifnot(exprs = {
	all.equal(cov2theta(S), theta)
	all.equal(theta2cov(theta), S)
})


## in_place_ragged_apply ###############################################

submean <- function(x) x - mean(x)

x <- 1:10
y <- in_place_ragged_apply(x, index = gl(2L, 5L), f = list(cumprod, submean))
f <- function(x) c(cumprod(x[1:5]), submean(x[6:10]))
stopifnot(all.equal(y, f(x)))

x <- as.data.frame(replicate(3L, c(exp(rnorm(5L)), qlogis(runif(5L)))))
y <- in_place_ragged_apply(x, index = gl(2L, 5L), f = list(log, plogis))
f <- function(x) c(log(x[1:5]), plogis(x[6:10]))
stopifnot(identical(y, as.data.frame(lapply(x, f))))
