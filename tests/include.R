library(epigrowthfit)
options(warn = 2L, error = recover)


tmp <- tempfile()
dir.create(tmp)
file.copy("src", tmp, recursive = TRUE)
setwd(file.path(tmp, "src"))

dll <- "test"
cpp <- paste0(dll, ".cpp")
TMB::compile(cpp)
dyn.load(TMB::dynlib(dll))

tt <- paste(readLines(cpp), collapse = "\n")
tt <- sub("^.*?enum[ \t\n]+test[ \t\n]*\\{(.*?)\\};.*$", "\\1", tt)
tt <- gsub("[ \t\n]", "", tt)
test_enums <- strsplit(tt, ",")[[1L]]
get_test_res <- function(test_enum, ...) {
    data <- list(test_flag = match(test_enum, test_enums, 0L) - 1L, ...)
    obj <- MakeADFun(data = data,
                     parameters = list(),
                     type = "Fun",
                     checkParameterOrder = FALSE,
                     DLL = dll)
    obj$report()$res
}

mvlgamma <- function(x, p) {
    0.25 * p * (p - 1) * log(pi) + rowSums(lgamma(outer(x, seq.int(0, p - 1, by = 1) / 2, `-`)))
}
dlkj <- function(x, eta, give_log = FALSE) {
    n <- 0.5 * (1 + sqrt(1 + 8 * length(x)))
    R <- diag(n)
    R[upper.tri(R)] <- x
    log_res <- (eta - 1) * (-sum(log(colSums(R * R))))
    if (give_log) log_res else exp(log_res)
}
dwishart <- function(x, df, scale, give_log = FALSE) {
    n <- 0.5 * (-1 + sqrt(1 + 8 * length(x)))
    X <- theta2cov(x)
    S <- theta2cov(scale)
    log_res <- -0.5 * (df * log(det(S)) + (-df + n + 1) * log(det(X)) + n * df * log(2) + 2 * mvlgamma(df / 2, n) + sum(diag(solve(S, X))))
    if (give_log) log_res else exp(log_res)
}
dinvwishart <- function(x, df, scale, give_log = FALSE) {
    n <- 0.5 * (-1 + sqrt(1 + 8 * length(x)))
    X <- theta2cov(x)
    S <- theta2cov(scale)
    log_res <- -0.5 * (-df * log(det(S)) + (df + n + 1) * log(det(X)) + n * df * log(2) + 2 * mvlgamma(df / 2, n) + sum(diag(solve(X, S))))
    if (give_log) log_res else exp(log_res)
}

test_that("list_of_vectors_t", {
    x <- list(rnorm(10L), seq_len(5L), TRUE, double(0L))
    res <- get_test_res(test_enum = "list_of_vectors_t", x = x)
    expect_identical(res, lapply(x, as.numeric))
})

test_that("is_NA_real_", {
    x <- c(0, NA, NaN, Inf)
    res <- get_test_res(test_enum = "is_NA_real_", x = x)
    expect_identical(res, c(0, 1, 0, 0))
})

test_that("is_finite", {
    x <- c(0, NA, NaN, Inf)
    res <- get_test_res(test_enum = "is_finite", x = x)
    expect_identical(res, c(1, 0, 0, 0))
})

test_that("logspace_diff", {
    log_x <- log(cumsum(rlnorm(10L)))
    res <- get_test_res(test_enum = "logspace_diff", log_x = log_x)
    expect_equal(res, log(diff(exp(log_x))))
})

test_that("mvlgamma", {
    x <- 1:12
    p1 <- 1L
    res1 <- get_test_res(test_enum = "mvlgamma", x = x, p = p1)
    expect_equal(res1, mvlgamma(x, p1))
    p2 <- 4L
    res2 <- get_test_res(test_enum = "mvlgamma", x = x, p = p2)
    expect_equal(res2, mvlgamma(x, p2))
})

test_that("dlkj", {
    n <- 4L
    x <- rnorm(n * (n - 1L) / 2)
    eta <- 2
    res <- get_test_res(test_enum = "dlkj", x = x, eta = eta, give_log = 1L)
    expect_equal(res, dlkj(x, eta, TRUE))
})

test_that("d(inv)?wishart", {
    n <- 4L
    x <- rnorm(n * (n + 1L) / 2)
    df <- 8
    scale <- rnorm(length(x))
    res1 <- get_test_res(test_enum = "dwishart",
                         x = x, df = df, scale = scale, give_log = 1L)
    expect_equal(res1, dwishart(x, df, scale, TRUE))
    res2 <- get_test_res(test_enum = "dinvwishart",
                         x = x, df = df, scale = scale, give_log = 1L)
    expect_equal(res2, dinvwishart(x, df, scale, TRUE))
})

test_that("dpois_robust", {
    log_lambda <- seq.int(0, 10, by = 1)
    x <- rpois(log_lambda, lambda = exp(log_lambda))
    res <- get_test_res(test_enum = "dpois_robust",
                        x = x, log_lambda = log_lambda, give_log = TRUE)
    expect_equal(res, dpois(x, lambda = exp(log_lambda), log = TRUE))
})

test_that("rnbinom_robust", {
    log_mu <- log(100)
    log_size <- log(50)
    n <- 1e+6L
    set.seed(10235L)
    res <- get_test_res(test_enum = "rnbinom_robust",
                        log_mu = log_mu, log_size = log_size, n = n)
    freq <- as.integer(table(factor(res, levels = seq.int(min(res), max(res)))))
    dens <- dnbinom(seq.int(min(res), max(res)),
                    mu = exp(log_mu), size = exp(log_size))
    expect_equal(freq / n, dens, tolerance = 1e-2)
})

test_that("eval_log_curve_exponential", {
    time <- seq.int(0, 100, by = 1)
    r <- log(2) / 20
    c0 <- 100

    res <- get_test_res(test_enum = "eval_log_curve_exponential",
                        time = time,
                        log_r = log(r),
                        log_c0 = log(c0))
    expect_equal(res, log(c0) + r * time)
})

test_that("eval_log_(curve|rt)_subexponential", {
    time <- seq.int(0, 100, by = 1)
    alpha <- log(2) / 20
    c0 <- 100
    p <- 0.95

    res1 <- get_test_res(test_enum = "eval_log_curve_subexponential",
                         time = time,
                         log_alpha = log(alpha),
                         log_c0 = log(c0),
                         logit_p = qlogis(p))
    expect_equal(res1, log(c0) + log1p((1 - p) * alpha * time / c0^(1 - p)) / (1 - p))

    res2 <- get_test_res(test_enum = "eval_log_rt_subexponential",
                         log_curve = res1,
                         log_alpha = log(alpha),
                         logit_p = qlogis(p))
    expect_equal(res2, log(alpha) - (1 - p) * res1)
})

test_that("eval_log_(curve|rt)_gompertz", {
    time <- seq.int(0, 100, by = 1)
    alpha <- log(2) / 20
    tinfl <- 100
    K <- 25000

    res1 <- get_test_res(test_enum = "eval_log_curve_gompertz",
                         time = time,
                         log_alpha = log(alpha),
                         log_tinfl = log(tinfl),
                         log_K = log(K))
    expect_equal(res1, log(K) - exp(-alpha * (time - tinfl)))

    res2 <- get_test_res(test_enum = "eval_log_rt_gompertz",
                         log_curve = res1,
                         log_alpha = log(alpha),
                         log_K = log(K))
    expect_equal(res2, log(alpha) + log(log(K) - res1))
})

test_that("eval_log_(curve|rt)_logistic", {
    time <- seq.int(0, 100, by = 1)
    r <- log(2) / 20
    tinfl <- 100
    K <- 25000

    res1 <- get_test_res(test_enum = "eval_log_curve_logistic",
                         time = time,
                         log_r = log(r),
                         log_tinfl = log(tinfl),
                         log_K = log(K))
    expect_equal(res1, log(K) - log1p(exp(-r * (time - tinfl))))

    res2 <- get_test_res(test_enum = "eval_log_rt_logistic",
                         log_curve = res1,
                         log_r = log(r),
                         log_K = log(K))
    expect_equal(res2, log(r) + log1p(-exp(res1) / K))
})

test_that("eval_log_(curve|rt)_richards", {
    time <- seq.int(0, 100, by = 1)
    r <- log(2) / 20
    tinfl <- 100
    K <- 25000
    a <- 1.005

    res1 <- get_test_res(test_enum = "eval_log_curve_richards",
                         time = time,
                         log_r = log(r),
                         log_tinfl = log(tinfl),
                         log_K = log(K),
                         log_a = log(a))
    expect_equal(res1, log(K) - log1p(a * exp(-a * r * (time - tinfl))) / a)

    res2 <- get_test_res(test_enum = "eval_log_rt_richards",
                         log_curve = res1,
                         log_r = log(r),
                         log_K = log(K),
                         log_a = log(a))
    expect_equal(res2, log(r) + log1p(-(exp(res1) / K)^a))
})

test_that("logspace_add_(baseline|offsets)", {
    time <- seq.int(0, 100, by = 1)
    log_curve <- log(100) + (log(2) / 20) * time
    log_diff_curve <- log(diff(exp(log_curve)))

    log_b <- log(2)
    res1 <- get_test_res(test_enum = "logspace_add_baseline",
                         log_curve = log_curve,
                         time = time,
                         log_b = log_b)
    expect_equal(res1, log(exp(log_b) * time + exp(log_curve)))

    log_w <- log(c(1, 0.9, 1.1, 0.8, 1.2, 0.7, 1.3))
    from <- 5L
    res2 <- get_test_res(test_enum = "logspace_add_offsets",
                         log_diff_curve = log_diff_curve,
                         log_w = log_w,
                         from = from)
    expect_equal(res2, log_diff_curve + rep_len(c(log_w[-seq_len(from)], log_w[seq_len(from)]), length(log_diff_curve)))
})

dyn.unload(TMB::dynlib(dll))
