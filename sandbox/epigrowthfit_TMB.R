library(TMB)
library(epigrowthfit)
library(tidyverse)
compile("epigrowthfit_TMB.cpp")
dyn.load("epigrowthfit_TMB.so")

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

## translate from x0 to t_half?
MakeADFun(data=list(t=dd$time, x=dd$x),
          parameters=with(i, list(x, ...)))
