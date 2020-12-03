library("TMB")
compile("test.cpp")
dyn.load(dynlib("test"))
set.seed(123)
m <- MakeADFun(data=list(x=1:10, y=1:10 + rnorm(10)),
               parameters=list(a=0, b=0, log_sigma=0),
               DLL="test",
               silent=TRUE)

r <- m$report(m$par)
r
