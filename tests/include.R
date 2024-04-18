library(epigrowthfit)
options(warn = 2L, error = if (interactive()) recover)

src <- "src"
src. <- system.file("tests", src, package = "epigrowthfit", mustWork = TRUE)
file.copy(c(src, src.), tempdir(), recursive = TRUE)
setwd(file.path(tempdir(), src))
Sys.setenv(R_TESTS = "") # startup.Rs does not exist here

dll <- "test"
cpp <- paste0(dll, ".cpp")
TMB::compile(cpp)
dyn.load(paste0(dll, .Platform[["dynlib.ext"]]))

getReport <- local({
	tt <- readLines(cpp)
	tt <- paste(tt, collapse = "\n")
	tt <- chartr("\t\r\n", "   ", tt)
	tt <- sub("^.*?enum +test *\\{(.*?)\\} *;.*$", "\\1", tt)
	tt <- gsub(" ", "", tt)
	enum. <- strsplit(tt, ",")[[1L]]

	function(enum, ...)
		TMB::MakeADFun(data = list(flag = match(enum, enum., 0L) - 1L, ...),
		               parameters = list(),
		               type = "Fun",
		               checkParameterOrder = FALSE,
		               DLL = dll)[["report"]]()[["ans"]]
})

mvlgamma <- function(x, p)
    0.25 * p * (p - 1) * log(pi) + rowSums(lgamma(outer(x, seq.int(0, by = 0.5, length.out = p), `-`)))

dlkj <- function(x, eta, log = FALSE) {
	n <- 0.5 * (1 + sqrt(1 + 8 * length(x)))
	R <- diag(n)
	R[upper.tri(R)] <- x
	log.ans <- (eta - 1) * (-sum(log(colSums(R * R))))
	if (log) log.ans else exp(log.ans)
}

dwishart <- function(x, df, scale, log = FALSE) {
	n <- 0.5 * (-1 + sqrt(1 + 8 * length(x)))
	X <- theta2cov(x)
	S <- theta2cov(scale)
	log.ans <- -0.5 * (df * log(det(S)) + (-df + n + 1) * log(det(X)) + n * df * log(2) + 2 * mvlgamma(0.5 * df, n) + sum(diag(solve(S, X))))
	if (log) log.ans else exp(log.ans)
}

dinvwishart <- function(x, df, scale, log = FALSE) {
	n <- 0.5 * (-1 + sqrt(1 + 8 * length(x)))
	X <- theta2cov(x)
	S <- theta2cov(scale)
	log.ans <- -0.5 * (-df * log(det(S)) + (df + n + 1) * log(det(X)) + n * df * log(2) + 2 * mvlgamma(0.5 * df, n) + sum(diag(solve(X, S))))
	if (log) log.ans else exp(log.ans)
}


## list_of_vectors_t
x <- list(rnorm(10L), seq_len(5L), TRUE, double(0L))
ans <- getReport("list_of_vectors_t", x = x)
stopifnot(identical(ans, lapply(x, as.double)))

## is_NA_real_
x <- c(0, NA, NaN, Inf)
ans <- getReport("is_NA_real_", x = x)
stopifnot(identical(ans, c(0, 1, 0, 0)))

## is_finite
x <- c(0, NA, NaN, Inf)
ans <- getReport("is_finite", x = x)
stopifnot(identical(ans, c(1, 0, 0, 0)))

## logspace_diff
log.x <- log(cumsum(rlnorm(10L)))
ans <- getReport("logspace_diff", log_x = log.x)
stopifnot(all.equal(ans, log(diff(exp(log.x)))))

## mvlgamma
x <- seq.int(1, 12, by = 1)
p1 <- 1L
ans1 <- getReport("mvlgamma", x = x, p = p1)
stopifnot(all.equal(ans1, mvlgamma(x, p1)))
p2 <- 4L
ans2 <- getReport("mvlgamma", x = x, p = p2)
stopifnot(all.equal(ans2, mvlgamma(x, p2)))

## dlkj
n <- 4L
x <- rnorm(0.5 * n * (n - 1L))
eta <- 2
ans <- getReport("dlkj", x = x, eta = eta, give_log = 1L)
stopifnot(all.equal(ans, dlkj(x, eta, TRUE)))

## d(inv)?wishart
n <- 4L
x <- rnorm(0.5 * n * (n + 1L))
df <- 8
scale <- rnorm(length(x))
ans1 <- getReport(   "dwishart", x = x, df = df, scale = scale, give_log = 1L)
stopifnot(all.equal(ans1,    dwishart(x, df, scale, TRUE)))
ans2 <- getReport("dinvwishart", x = x, df = df, scale = scale, give_log = 1L)
stopifnot(all.equal(ans2, dinvwishart(x, df, scale, TRUE)))

## dpois_robust
log.lambda <- seq.int(0, 10, by = 1)
x <- rpois(log.lambda, lambda = exp(log.lambda))
ans <- getReport("dpois_robust", x = x, log_lambda = log.lambda, give_log = TRUE)
stopifnot(all.equal(ans, dpois(x, lambda = exp(log.lambda), log = TRUE)))

## rnbinom_robust
log.mu <- log(100)
log.size <- log(50)
n <- 1e+06L
set.seed(10235L)
ans <- getReport("rnbinom_robust", log_mu = log.mu, log_size = log.size, n = n)
mm <- seq.int(min(ans), max(ans))
freq <- as.integer(table(factor(ans, levels = mm)))
dens <- dnbinom(mm, mu = exp(log.mu), size = exp(log.size))
stopifnot(all.equal(freq / n, dens, tolerance = 1e-02))

## eval_log_curve_exponential
time <- seq.int(0, 100, by = 1)
r <- log(2) / 20
c0 <- 100
ans <- getReport("eval_log_curve_exponential",
                 time = time,
                 log_r = log(r),
                 log_c0 = log(c0))
stopifnot(all.equal(ans, log(c0) + r * time))

## eval_log_(curve|rt)_subexponential
time <- seq.int(0, 100, by = 1)
alpha <- log(2) / 20
c0 <- 100
p <- 0.95
ans1 <- getReport("eval_log_curve_subexponential",
                  time = time,
                  log_alpha = log(alpha),
                  log_c0 = log(c0),
                  logit_p = qlogis(p))
all.equal(ans1, log(c0) + log1p((1 - p) * alpha * time / c0^(1 - p)) / (1 - p))
ans2 <- getReport("eval_log_rt_subexponential",
                  log_curve = ans1,
                  log_alpha = log(alpha),
                  logit_p = qlogis(p))
stopifnot(all.equal(ans2, log(alpha) - (1 - p) * ans1))

## eval_log_(curve|rt)_gompertz
time <- seq.int(0, 100, by = 1)
alpha <- log(2) / 20
tinfl <- 100
K <- 25000
ans1 <- getReport("eval_log_curve_gompertz",
                  time = time,
                  log_alpha = log(alpha),
                  log_tinfl = log(tinfl),
                  log_K = log(K))
stopifnot(all.equal(ans1, log(K) - exp(-alpha * (time - tinfl))))
ans2 <- getReport("eval_log_rt_gompertz",
                  log_curve = ans1,
                  log_alpha = log(alpha),
                  log_K = log(K))
stopifnot(all.equal(ans2, log(alpha) + log(log(K) - ans1)))

## eval_log_(curve|rt)_logistic
time <- seq.int(0, 100, by = 1)
r <- log(2) / 20
tinfl <- 100
K <- 25000
ans1 <- getReport("eval_log_curve_logistic",
                  time = time,
                  log_r = log(r),
                  log_tinfl = log(tinfl),
                  log_K = log(K))
stopifnot(all.equal(ans1, log(K) - log1p(exp(-r * (time - tinfl)))))
ans2 <- getReport("eval_log_rt_logistic",
                  log_curve = ans1,
                  log_r = log(r),
                  log_K = log(K))
stopifnot(all.equal(ans2, log(r) + log1p(-exp(ans1) / K)))

## eval_log_(curve|rt)_richards
time <- seq.int(0, 100, by = 1)
r <- log(2) / 20
tinfl <- 100
K <- 25000
a <- 1.005
ans1 <- getReport("eval_log_curve_richards",
                  time = time,
                  log_r = log(r),
                  log_tinfl = log(tinfl),
                  log_K = log(K),
                  log_a = log(a))
stopifnot(all.equal(ans1, log(K) - log1p(a * exp(-a * r * (time - tinfl))) / a))
ans2 <- getReport("eval_log_rt_richards",
                  log_curve = ans1,
                  log_r = log(r),
                  log_K = log(K),
                  log_a = log(a))
stopifnot(all.equal(ans2, log(r) + log1p(-(exp(ans1) / K)^a)))

## logspace_add_(baseline|offsets)
time <- seq.int(0, 100, by = 1)
log_curve <- log(100) + (log(2) / 20) * time
log_diff_curve <- log(diff(exp(log_curve)))
log_b <- log(2)
log_w <- log(c(1, 0.9, 1.1, 0.8, 1.2, 0.7, 1.3))
from <- 5L
ans1 <- getReport("logspace_add_baseline",
                  log_curve = log_curve,
                  time = time,
                  log_b = log_b)
stopifnot(all.equal(ans1, log(exp(log_b) * time + exp(log_curve))))
ans2 <- getReport("logspace_add_offsets",
                  log_diff_curve = log_diff_curve,
                  log_w = log_w,
                  from = from)
stopifnot(all.equal(ans2, log_diff_curve + rep_len(c(log_w[-seq_len(from)], log_w[seq_len(from)]), length(log_diff_curve))))
