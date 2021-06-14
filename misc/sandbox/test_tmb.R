library("TMB")
compile("test_tmb.cpp")
dyn.load(dynlib("test_tmb"))
set.seed(123)
m <- MakeADFun(data=list(x=1:10, y=1:10 + rnorm(10)),
               parameters=list(a=0, b=0, log_sigma=0),
               DLL="test_tmb",
               silent=TRUE)

r <- m$report(m$par)
r
