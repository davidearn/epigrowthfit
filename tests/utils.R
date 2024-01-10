library(epigrowthfit)
options(warn = 2L, error = recover)


test_that("disambiguate", {
    x <- c("c", "b", "a", "b", "a",
           "a", "b", "c", "c", "c")
    y <- c("c[1]", "b[1]", "a[1]", "b[2]", "a[2]",
           "a[3]", "b[3]", "c[2]", "c[3]", "c[4]")
    expect_identical(disambiguate(x), y)

    x <- `names<-`(seq_along(x), x)
    y <- `names<-`(seq_along(x), y)
    expect_identical(disambiguate(x, nms = TRUE), y)
})

test_that("literal_rle", {
    x <- c(0, NA, NaN, Inf, 1)
    times <- 1:5
    y <- rep.int(x, times)
    rle_y <- literal_rle(y)
    expect_identical(rle_y, list(lengths = times, values = x))
    expect_identical(y, inverse.rle(rle_y))
})

test_that("locf", {
    x <- c(NA, NA, 1, NA, 2, 2, 3, NA)
    expect_identical(locf(x), c(NA, NA, 1, 1, 2, 2, 3, 3))
    expect_identical(locf(x, x0 = 0), c(0, 0, 1, 1, 2, 2, 3, 3))
})

test_that("wald", {
    estimate <- rnorm(6L, 0, 1)
    se <- rlnorm(6L, 0, 0.1)
    level <- 0.95
    q <- qchisq(level, df = 1)
    W <- wald(estimate = estimate, se = se, level = level)
    expect_type(W, "double")
    expect_true(is.matrix(W))
    expect_identical(dim(W), c(6L, 2L))
    expect_identical(dimnames(W), list(NULL, c("lower", "upper")))
    expect_equal(W[, "lower"], estimate - sqrt(q) * se)
    expect_equal(W[, "upper"], estimate + sqrt(q) * se)
})

test_that("cov2theta, theta2cov", {
    set.seed(230719L)
    n <- 5L
    S <- crossprod(matrix(rnorm(n * n), n, n))
    R <- chol(S)
    R <- R * rep.int(1 / diag(R), rep.int(n, n))
    theta <- c(0.5 * log(diag(S)), R[upper.tri(R)])
    expect_equal(cov2theta(S), theta)
    expect_equal(theta2cov(theta), S)
})

test_that("in_place_ragged_apply", {
    submean <- function(x) x - mean(x)

    x <- 1:10
    y <- in_place_ragged_apply(x, index = gl(2L, 5L), f=list(cumprod, submean))
    f <- function(x) c(cumprod(x[1:5]), submean(x[6:10]))
    expect_equal(y, f(x))

    x <- as.data.frame(replicate(3L, c(exp(rnorm(5L)), qlogis(runif(5L)))))
    y <- in_place_ragged_apply(x, index = gl(2L, 5L), f = list(log, plogis))
    f <- function(x) c(log(x[1:5]), plogis(x[6:10]))
    expect_identical(y, as.data.frame(lapply(x, f)))
})
