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


if (FALSE) {
## rle1 ################################################################

x <- c(0, NA_real_, NaN, Inf, 1)
times <- 1:5
y <- rep.int(x, times)
rle.y <- rle1(y)
stopifnot(exprs = {
	identical(rle.y, list(lengths = times, values = x))
	identical(y, inverse.rle(rle.y))
})


## locf ################################################################

x <- c(NA, NA, 1, NA, 2, 2, 3, NA)
stopifnot(exprs = {
	identical(locf(x), c(NA, NA, 1, 1, 2, 2, 3, 3))
	identical(locf(x, x0 = 0), c(0, 0, 1, 1, 2, 2, 3, 3))
})
}


## wald ################################################################

value <- rnorm(6L, 0, 1)
se <- rlnorm(6L, 0, 0.1)
level <- 0.95
q <- qchisq(level, df = 1)
W <- wald(value = value, se = se, level = level)
stopifnot(exprs = {
	is.double(W)
	identical(dim(W), c(6L, 2L))
	is.null(dimnames(W))
	all.equal(W[, 1L], value - sqrt(q) * se)
	all.equal(W[, 2L], value + sqrt(q) * se)
})


## cov2theta, theta2cov ################################################

set.seed(230719L)
n <- 6L
Sigma <- crossprod(matrix(rnorm(n * n), n, n))
R <- chol(Sigma)
R1 <- R * rep(1 / diag(R, names = FALSE), each = n)
theta <- c(0.5 * log(diag(Sigma, names = FALSE)), R1[upper.tri(R1)])
stopifnot(exprs = {
	all.equal(cov2theta(Sigma), theta)
	all.equal(theta2cov(theta), Sigma)
})
