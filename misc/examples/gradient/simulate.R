library("epigrowthfit")

model <- egf_model(curve = "logistic", family = "nbinom")

r <- log(2) / 20
tinfl <- 160
K <- 25000
nbdisp <- 50

mu <- log(c(r, tinfl, K, nbdisp))
Sigma <- diag(rep_len(1, length(mu)))

zz <- simulate(model,
  nsim = 100L,
  seed = 202737L,
  mu = mu,
  Sigma = Sigma,
  cstart = 10
)
mm0 <- egf(zz, do_fit = FALSE)

Y <- mm0$Y_init
beta <- as.numeric(colMeans(Y))
b <- as.numeric(t(Y) - beta)
theta <- rep_len(0, 10L)
theta[1:3] <- log(apply(Y[, 1:3], 2L, sd))

l <- list(
  actual = zz$actual,
  data = mm0$tmb_args$data[c("time", "time_seg_len", "x", "Z")],
  parameters = list(beta = beta, b = b, theta = theta)
)
saveRDS(l, file = "gradient.rds")
