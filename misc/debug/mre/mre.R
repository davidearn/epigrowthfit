library("TMB")
packageVersion("TMB") # 1.7.20

zz <- "mre"
compile(paste0(zz, ".cpp"))
dyn.load(dynlib(zz))
load(paste0(zz, ".RData")) # actual, tmb_data, tmb_parameters

openmp(4L)
obj <- MakeADFun(
  data = tmb_data,
  parameters = tmb_parameters,
  random = "b",
  DLL = zz,
  inner.method = "newton",
  inner.control = list(maxit = 1000L, trace = 1L),
  silent = FALSE
)
opt <- with(obj, optim(par, fn, gr, method = "BFGS", control = list(maxit = 1000L, trace = 1L)))

zz <- data.frame(
  ## Generative model (used to simulate data)
  actual,
  ## Fitted model
  fitted = obj$env$last.par.best
)
zz[-obj$env$random, ]

## Maximum gradient component
max(abs(obj$gr(obj$env$last.par.best[-obj$env$random])))
