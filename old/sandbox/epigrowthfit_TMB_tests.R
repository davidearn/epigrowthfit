## exponential growth/Poisson
library(tidyverse)
library(TMB)

compile("epigrowthfit_TMB.cpp")
dyn.load(dynlib("epigrowthfit_TMB"))

### 

set.seed(101)
dd <- (tibble(time=1:30,
             cum_inc=10*exp(0.1*time),
             inc=c(NA,diff(cum_inc)),
             x=rpois(30,lambda=inc)
             )
    %>% drop_na(x)
)

## FIXME: make sure epigrowthfit is robust against NAs?

plot(x~time,data=dd)

eg_TMB <- function(data,distrib="dbinom",model="richards") {
    ## data <- dd
    i <- with(data,
              getInits(time, x, model="richards", distrib="nbinom"))
    dd <- data[i$first:i$peak,]
    all_params <- with(i, list(log_thalf=log(tail(dd$time,1))
                         , log_x0 = log(theta0[["x0"]])
                         , log_K=log(theta0[["K"]])
                         , log_r=log(theta0[["r"]])
                         , log_p=log(theta0[["s"]])
                         , log_nbdisp=log(theta0[["ll.k"]])
                           ))
    dt <- diff(dd$time[1:2]) ## assume equal time steps
    ad_data <- list(t=c(dd$time[1]-dt,dd$time), x=dd$x, debug=0,
                    curve_flag=3, distr_flag=2)
    m <- MakeADFun(data=ad_data, parameters=all_params,
                   DLL="epigrowthfit_TMB",
                   map=list(log_x0=factor(NA)),
                   silent=TRUE)
    params <- all_params[!names(all_params) %in% "log_x0"]
    fit <- nlminb(start=unlist(params), objective=m$fn, gradient=m$gr)
    fit$objective
    sdreport(m,m$env$last.par.best)
    ad_data_pois <- ad_data
    ad_data_pois$distr_flag <- 1
    m_pois <- MakeADFun(data=ad_data_pois, parameters=all_params,
                   DLL="epigrowthfit_TMB",
                   map=list(log_x0=factor(NA), log_nbdisp=factor(NA)),
                   silent=TRUE)
    params_pois <- all_params[!names(all_params) %in% c("log_x0","log_nbdisp")]
    fit_pois <- nlminb(start=unlist(params_pois), objective=m_pois$fn, gradient=m_pois$gr)
    sdreport(m_pois,m_pois$env$last.par.best)
}
    


