data("sleepstudy",package="lme4")
library(glmmTMB)
library(tmbstan)
fm1 <- glmmTMB(Reaction ~ Days + (Days|Subject),
               sleepstudy)
m <- tmbstan(fm1$obj,cores=4)
dd <- as.data.frame(m)
str(dd)
compval <- with(dd,`b[1]`*exp(`theta[1]`)+`beta[1]`)
quantile(compval, c(0.025,0.975))

library(epigrowthfit)
library(rstan)
example(egf)
m2 <- tmbstan(y1$madf_out,cores=4)
names(dd2 <- as.data.frame(m2))
## only chain 1 makes any sense at all ... these are first 1000 samples ...
cbind(coef(y1,log=TRUE),
      colMeans(dd2[1:1000,1:4]))
traceplot(m2)
library(ggplot2)
library(tidyr)
library(dplyr)
dd2w <- dd2 %>% mutate(n=seq(nrow(dd2))) %>%
    pivot_longer(cols=-n) %>% filter(name !="lp__", n<=1000)
ggplot(dd2w,aes(n,value)) + geom_line() + facet_wrap(~name, scale="free")
ggplot(dd2w,aes(value)) + geom_histogram(bins=50) + facet_wrap(~name, scale="free")

## don't use crazy starting value
## pick starting conditions from a uniform cube of (by default) +/- 30% around
##  true coef values
ifun <- function(m,s=0.3) {
    function() {
        cc <- coef(m, log=TRUE)
        return(as.list(cc * runif(length(cc), 1-s, 1+s)))
    }
}

set.seed(101)
m3 <- tmbstan(y1$madf_out,
              init=ifun(y1),
              cores=4)
traceplot(m3)
