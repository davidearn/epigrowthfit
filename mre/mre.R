library("TMB")
compile("mre.cpp")
dyn.load(dynlib("mre"))
load("mre.RData")

# outfile <- file("mre.Rout", open = "wt")
# sink(outfile, type = "output")
# sink(outfile, type = "message")

obj <- MakeADFun(
  data = tmb_data,
  parameters = tmb_parameters,
  random = "b",
  DLL = "mre",
  inner.method = "newton",
  inner.control = list(maxit = 1000L),
  silent = FALSE
)
opt <- with(obj, nlminb(par, fn, gr, control = list(trace = 1L)))
sdr <- try(sdreport(obj))

# sink(type = "message")
# sink(type = "output")
