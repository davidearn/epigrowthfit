library("TMB")
packageVersion("TMB") # 1.7.20
dll <- "gradient"

l <- readRDS(paste0(dll, ".rds")) # list(data, parameters)
compile(paste0(dll, ".cpp"), openmp = TRUE)
dyn.load(dynlib(dll))

openmp(4L)
obj <- MakeADFun(
  data = l$data,
  parameters = l$parameters,
  random = "b",
  DLL = dll,
  inner.control = list(maxit = 1000L)
)

p0 <- obj$par
n <- length(p0)
res <- matrix(NA_real_, nrow = n + 1L + n, ncol = 3L)

set.seed(110108L)
for (i in seq_len(ncol(res))) {
  p1 <- optim(p0, obj$fn, obj$gr, method = "BFGS", control = list(maxit = 1000L, trace = 1L))$par
  f1 <- obj$fn(p1)
  g1 <- obj$gr(p1)
  if (length(g1) != n) {
    g1 <- rep_len(NaN, n)
  }
  res[, i] <- c(p1, f1, g1)
  dp <- abs(p1 - p0)
  p0 <- p1 + rnorm(n, 0, pmin(1e-03, pmax(1e-06, 0.25 * dp)))
}
r <- sdreport(obj)
