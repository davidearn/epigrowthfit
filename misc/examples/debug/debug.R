library("epigrowthfit")
options(warn = 1L, egf.cores = 4L)

## Simulate, estimate, return `c(error, gradient, convergence)`
f <- function(model, N, mu, sigma) {
  zz <- simulate(model,
    nsim = N,
    mu = mu,
    Sigma = diag(rep_len(sigma^2, 4L)),
    cstart = 10
  )
  mm0 <- egf(zz, do_fit = FALSE)

  p <- length(mu)
  beta_init <- colMeans(mm0$Y_init)
  b_init <- t(mm0$Y_init) - beta_init
  theta_init <- rep_len(0, p * (p + 1) / 2)
  theta_init[seq_len(p - 1L)] <- log(apply(b_init[-p, ], 1L, sd))
  b_init <- b_init / exp(theta_init[seq_len(p)])

  mm <- try(update(mm0,
    do_fit = TRUE,
    control = egf_control(
      optimizer = egf_optimizer(f = optim, args = list(method = "BFGS"), control = list(maxit = 1000L, trace = FALSE)),
      inner_optimizer = egf_inner_optimizer(args = list(maxit = 1000L, trace = FALSE)),
      trace = FALSE
    ),
    init = c(beta_init, b_init, theta_init)
  ))

  if (inherits(mm, "try-error")) {
    cat("Failed to fit\n")
    return(rep_len(NaN, p * (p + 3L) + 1L))
  }
  v1 <- (mm$best - zz$actual)[mm$nonrandom]
  v2 <- mm$tmb_out$gr(mm$best[mm$nonrandom])
  v3 <- mm$optimizer_out$convergence
  if (!(is.numeric(v2) && length(v2) == length(v1) && all(is.finite(v2)))) {
    cat("Failed to compute outer gradient\n")
    v2 <- v1
    v2[] <- NaN
  }
  if (v3 != 0) {
    cat("Nonzero convergence status\n")
  }
  unname(c(v1, v2, v3))
}

## Model
model <- egf_model(curve = "logistic", family = "nbinom")

## Number of time series per simulation
N <- 300L

## Number of simulations
n <- 30L

## Nonlinear and dispersion model parameters
r <- log(2) / 20
tinfl <- 160
K <- 25000
nbdisp <- 50

## Means
mu <- log(c(r, tinfl, K, nbdisp))
p <- length(mu)

## Standard deviations across time series
sigma <- c(0.01, 0.05, 0.1, 0.5, 1)

res <- array(NA_real_, dim = c(p * (p + 3L) + 1L, n, length(sigma)))
set.seed(154051L)
for (k in seq_along(sigma)) {
  for (j in seq_len(n)) {
    cat("commencing simulation", j, "of", n, "for sd value", k, "of", length(sigma), "\n")
    res[, j, k] <- f(model = model, N = N, mu = mu, sigma = sigma[k])
  }
}
dn1 <- sprintf("%s_%s[%d]",
  rep(c("error", "gradient"), each = p * (p + 3) / 2),
  rep(c("beta", "theta"), times = c(p, p * (p + 1) / 2)),
  c(seq_len(p), seq_len(p * (p + 1) / 2))
)
dimnames(res) <- list(
  c(dn1, "convergence"),
  NULL,
  sprintf("sd=%.2f", sigma)
)
save.image("debug1.RData")


# load("debug.RData")
# pdf("stripchart_errors.pdf", width = 8, height = 6)
#
# par(mfcol = c(3L, 5L), mar = c(4.5, 4, 1, 1) + 0.1,
#     mgp = c(3, 0.7, 0), las = 2, cex.main = 0.9)
# for (k in seq_along(sigma)) {
#   do_title <- function() {
#     title(main = sprintf("sigma = %.2f", sigma[k]), adj = 0, line = 0.25)
#   }
#
#   dd <- data.frame(
#     x = as.numeric(res[seq_len(p), , k]),
#     g = gl(p, 1L, p * n, sprintf("beta[%d]", seq_len(p))),
#     stringsAsFactors = TRUE
#   )
#   stripchart(x ~ g, data = dd, vertical = TRUE, ylab = "error",
#              ylim = c(-0.4, 0.4), pch = 16, col = "#00000054")
#   do_title()
#
#   dd <- data.frame(
#     x = as.numeric(res[p + seq_len(p), , k]),
#     g = gl(p, 1L, p * n, sprintf("theta[%d]", p)),
#     stringsAsFactors = TRUE
#   )
#   stripchart(x ~ g, data = dd, vertical = TRUE, ylab = "error",
#              ylim = c(-15, 5), pch = 16, col = "#00000054")
#   do_title()
#
#   l <- p * (p - 1) / 2
#   dd <- data.frame(
#     x = as.numeric(res[2 * p + seq_len(l), , k]),
#     g = gl(l, 1L, l * n, sprintf("theta[%d]", 2 * p + seq_len(l))),
#     stringsAsFactors = TRUE
#   )
#   stripchart(x ~ g, data = dd, vertical = TRUE, ylab = "error",
#              ylim = c(-20, 20), pch = 16, col = "#00000054")
#   do_title()
# }
# dev.off()
#
# zz <- apply(res[q + seq_len(q), , ], 2:3, function(x) max(abs(x)))
# zz <- apply(zz, 2L, sort, decreasing = TRUE)
# colnames(zz) <- sprintf("sigma=%.2f", sigma)
