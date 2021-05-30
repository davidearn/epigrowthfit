# https://github.com/kaskr/adcomp/issues/348

library("TMB")
# packageVersion("TMB") # 1.7.20

zz <- "mre"
compile(paste0(zz, ".cpp"))
dyn.load(dynlib(zz))
load(paste0(zz, ".RData")) # tmb_data, tmb_parameters

# outfile <- file(paste0(zz, ".Rout"), open = "wt")
# sink(outfile, type = "output")
# sink(outfile, type = "message")

obj <- MakeADFun(
  data = tmb_data,
  parameters = tmb_parameters,
  random = "b",
  DLL = zz,
  inner.method = "newton",
  inner.control = list(maxit = 1000L),
  silent = FALSE
)
opt <- with(obj, nlminb(par, fn, gr, control = list(trace = 1L)))
sdr <- try(sdreport(obj))

# sink(type = "message")
# sink(type = "output")
