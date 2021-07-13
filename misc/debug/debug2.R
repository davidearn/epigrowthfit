library("epigrowthfit")

N <- 300L
model <- egf_model(curve = "logistic", family = "nbinom")

r <- log(2) / 20
tinfl <- 100
K <- 25000
nbdisp <- 50

mu <- log(c(r, tinfl, K, nbdisp))
p <- length(mu)
Sigma <- diag(rep_len(1, p))

set.seed(101958L)
zz <- egf_simulate(N = N, model = model, mu = mu, Sigma = Sigma, cstart = 10)
mm0 <- egf(zz, do_fit = FALSE)

beta_init <- colMeans(mm0$Y_init)
b_init <- t(mm0$Y_init) - beta_init
theta_init <- c(log(apply(b_init, 1L, sd)), rep_len(0, 6L))
theta_init[p] <- 0 # replacing log(0) = -Inf
b_init <- b_init / exp(theta_init[seq_len(p)])

mm <- update(mm0,
  do_fit = TRUE,
  control = egf_control(
    optimizer = egf_optimizer(f = optim, args = list(method = "BFGS"), control = list(maxit = 1000L, trace = FALSE)),
    inner_optimizer = egf_inner_optimizer(args = list(maxit = 1000L, trace = FALSE)),
    trace = FALSE,
    omp_num_threads = 4L
  ),
  init = c(beta_init, b_init, theta_init)
)

data.frame(actual = zz$actual, fitted = mm$best)[mm$nonrandom, ]
max(abs(mm$tmb_out$gr(mm$best[mm$nonrandom])))
