library("TMB")
packageVersion("TMB") # 1.7.20

rds <- list.files(pattern = "\\.rds$")
stopifnot(length(rds) == 1L)
l <- readRDS(rds) # list(actual, data, parameters)

cpp <- list.files(pattern = "\\.cpp$")
stopifnot(length(cpp) == 1L)
compile(cpp)

dll <- sub("\\.cpp$", "", cpp)
dyn.load(dynlib(dll))

## openmp(6L)
obj <- MakeADFun(
  data = l$data,
  parameters = l$parameters,
  random = "b",
  DLL = dll,
  inner.method = "newton",
  inner.control = list(maxit = 1000L, trace = 1L),
  silent = FALSE
)
opt <- with(obj, optim(par, fn, gr, method = "BFGS", control = list(maxit = 1000L, trace = 1L)))
opt$convergence # O

## Optimizer seems to have arrived somewhere reasonable
d <- data.frame(
  ## Generative model from which data were simulated
  actual = l$actual,
  ## Fitted model
  fitted = obj$env$last.par.best
)
random <- obj$env$random
d[-random, ]

## Yet the objective function gradient is extremely large
par <- d$fitted[-random]
(f0 <- obj$fn(par))
(g0 <- obj$gr(par))
max(abs(g0)) # roughly 2e+12

## Sampling a neighbourhood of `par` ...
set.seed(230902L)
res <- replicate(100L, {
  p <- par + rnorm(par, 0, 0.001)
  c(obj$fn(p), obj$gr(p))
})
all(res[1L, ] > f0) # TRUE
max(abs(res[-1L, ])) # roughly 1e+05
