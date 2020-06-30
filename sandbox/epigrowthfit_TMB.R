library(TMB)
library(epigrowthfit)
library(tidyverse)
library(bbmle)
compile("epigrowthfit_TMB.cpp")
dyn.load(dynlib("epigrowthfit_TMB"))

ontario <- (dplyr::filter(canada_covid,Province=="ON")
    ## need time to be numeric
    %>% mutate(time=lubridate::decimal_date(Date))
    %>% select(time, x=newConfirmations)
    ## NAs mess us up
    %>% drop_na(x)
)

i <- with(ontario,
          getInits(time, x, model="richards", distrib="nbinom"))

dd <- (ontario[i$first:i$peak,]
    %>% mutate_at("time", ~ . -min(.))
)

params <- with(i, list(log_thalf=log(tail(dd$time,1))
                     , log_K=log(theta0[["K"]])
                     , log_r=log(theta0[["r"]])
                     , log_p=log(theta0[["s"]])
                     , log_nb_disp=log(theta0[["ll.k"]])
))
## translate from x0 to t_half?
m <- MakeADFun(data=list(t=dd$time, x=dd$x, debug=0),
          parameters=params,
          DLL="epigrowthfit_TMB",
          silent=TRUE)
parnames(m$fn) <- names(params)
m$fn(unlist(params))
m$gr(unlist(params))

nlminb(start=unlist(params), objective=m$fn, gradient=m$gr)
mle2(minuslogl=m$fn, start=params, gr=m$gr)

plot(x~time,data=dd,log="y")
lines(dd$time,m$report()$inccurve)
lines(dd$time,m$report(unlist(params))$inccurve,col=2)
with(sdreport(m),matlines(dd$time,cbind(value,value-1.96*sd,
                                        value+1.96*sd),lty=c(1,2,2),
                          col=1))

## FIXME: first and last point are wonky because of differencing etc.?
