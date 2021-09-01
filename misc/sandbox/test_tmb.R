library("TMB")
packageVersion("TMB")
R.version.string

dll <- "test_tmb"
compile(paste0(dll, ".cpp"))
dyn.load(dynlib(dll))

res <- MakeADFun(
  data = list(),
  parameters = list(),
  type = "Fun",
  checkParameterOrder = FALSE,
  DLL = dll
)$report()
