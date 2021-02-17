
library(epigrowthfit)
library(dplyr)
library(lubridate)

## covid growth in Ontario
## https://wzmli.github.io/COVID19-Canada/COVID19_Canada.csv
canada.raw <- read.csv("canada-covid/COVID19_Canada.csv")

onc <- (canada.raw
    %>% filter(Province == "ON")
    %>% select(-Province)
    %>% mutate(date = as.Date(Date))
    %>% mutate(time = decimal_date(date))
    %>% mutate(cases = c(0,diff(confirmed_positive)))
    ## epigrowthfit does not handle NAs:
    %>% mutate(cases = na.omit(cases)) ## FIX: is this doing anything?
    %>% mutate(deaths = c(0,diff(deceased)))
    %>% select(date, time, cases, Hospitalization, deaths)
)

## naive fit:
onc.init <- egf_init(date = onc$date, cases = onc$cases,
                     curve="logistic")
onc.fit <- egf(onc.init)
## cumulative incidence
plot(onc.init)
plot(onc.fit)

## deaths
onc.deaths <- subset(onc, !is.na(deaths) & deaths >= 0)
onc.deaths.init <- egf_init(
    date = onc.deaths$date, cases = onc.deaths$deaths,
    curve="logistic")
onc.deaths.fit <- egf(onc.deaths.init)
plot(onc.deaths.init, inc="interval")
plot(onc.deaths.fit, inc="interval")

## hospitalization data (occupancy, not admissions)
## doesn't make sense to fit to occupancy
onc.hosp <- subset(onc, !is.na(Hospitalization))
onc.hosp.init <- egf_init(
    date = onc.hosp$date, cases = onc.hosp$Hospitalization,
    curve="logistic")
onc.hosp.fit <- egf(onc.hosp.init)
plot(onc.hosp.init, inc="interval")
plot(onc.hosp.fit, inc="interval", ylim=c(10,2000))

## first wave:
onc.init1 <- update(onc.init, peak=55)
onc.fit1 <- egf(onc.init1)
## cumulative incidence
plot(onc.init1)
plot(onc.fit1)
## interval incidence
try(plot(onc.init1, inc="interval"))
try(plot(onc.fit1, inc="interval"))

## second wave:
onc.init2 <- update(onc.init, first=205, peak=240)
onc.fit2 <- egf(onc.init2)
## cumulative incidence
plot(onc.init2)
plot(onc.fit2)
## interval incidence
try(plot(onc.init2, inc="interval"))
try(plot(onc.fit2, inc="interval"))

## second wave, part 2:
onc.init3 <- update(onc.init
                  , curve="exponential"
                  , first=255)
onc.fit3 <- egf(onc.init3
              ##, method="Nelder-Mead"
                )
## cumulative incidence
plot(onc.init3)
plot(onc.fit3)
## interval incidence
try(plot(onc.init3, inc="interval"))
try(plot(onc.fit3, inc="interval"))

show_doubling_time <- function(fit, h=0.5, pos=4, digits=2,
                               xshift = 0, ...) {
    xpos <- fit$init$last + xshift
    ycum <- fit$eval_cum_inc(fit$init$time)
    y <- c(0,diff(ycum))
    ymin <- y[fit$init$first]
    ymax <- y[fit$init$last]
    ypos <- ymin + h*(ymax-ymin)
    doubtime <- signif(compute_doubling_time(fit), digits)
    text(x=xpos, y=ypos, pos=pos, xpd=NA, ...,
         labels=sprintf("doubling time:\n%g days",
                 doubtime))
}

get_date <- function(x) {
    max(x$data$date)
}

show_date  <- function(fit, where="topleft", cex=0.8, ...) {
    legend(where, bty="n", cex=cex,
           legend=as.character(get_date(fit)))
}

## BOTH WAVES

## plot data only
if (!interactive()) pdf("ontario_covid.pdf", height=5)
plot(onc.fit1, inc="interval"
   , annotate=FALSE
   , ylab="Daily reported cases"
   , line_style=list(col=0)
   , polygon_style=list(col=0)
   ##, main="Ontario COVID-19 confirmed cases"
   , main=paste0("Ontario COVID-19 confirmed cases up to ", get_date(onc.fit1))
     )
if (!interactive()) dev.off()

ontario_fits <- function(log=TRUE, text_dbl3=list(x=245,y=500),
            style.pt=list(points_main=list(bg="darkred",col="black"))
                         ) {
## plot with fits
plot(onc.fit1, inc="interval"
   , annotate=FALSE
   , log=log
   , ylab="Daily reported cases"
   , style=style.pt
   ##, polygon_style=list(col="yellow", border="black")
   ##, main="Ontario COVID-19 confirmed cases"
   , main=paste0("Ontario COVID-19 confirmed cases up to ", get_date(onc.fit1))
     )
##show_doubling_time(onc.fit1, h=0.1, pos=4, xshift=-8)
plot(onc.fit2, inc="interval"
   , style=style.pt
   ##, polygon_style=list(col="yellow", border="black")
   , add=TRUE)
##show_doubling_time(onc.fit2, pos=4, h=0, xshift=-15)

## show_date(onc.fit2)
plot(onc.fit3, inc="interval"
   ##, polygon_style=list(col="yellow", border="black")
   , style=list(text_dbl=text_dbl3,
                    points_main=style.pt[["points_main"]])
   , add=TRUE)
}

if (!interactive()) pdf("ontario_covid_linear_fit.pdf", height=5)
ontario_fits(log=FALSE)
if (!interactive()) dev.off()
if (!interactive()) pdf("ontario_covid_growth_fit.pdf", height=5)
ontario_fits()
if (!interactive()) dev.off()

##################################################################

## WORLD COVID

## cases
owid.raw <- read.csv("owid-covid-daily-counts/new_cases.csv")
owid <- (owid.raw
    %>% mutate(date = as.Date(date)) 
    ##%>% mutate(cases = c(0,diff(World)))
    %>% rename(cases = World)
    %>% mutate(time = decimal_date(date))
    %>% select(date, time, cases)
)

owid.init <- egf_init(date = owid$date, cases = owid$cases,
                     curve="logistic")
owid.fit <- egf(owid.init)
## cumulative incidence
plot(owid.init)
plot(owid.fit)

if (!interactive()) pdf("world_covid.pdf", height=5)
plot(owid.fit, inc="interval"
   , ylim=c(10^2,10^6)
   , annotate=FALSE
   , ylab="Daily reported cases"
   , line_style=list(col=0)
   , polygon_style=list(col=0)
     ##, main="Worldwide COVID-19 confirmed cases"
   , main=paste0("Worldwide COVID-19 confirmed cases up to ", get_date(owid.fit))
)
if (!interactive()) dev.off()
if (!interactive()) pdf("world_covid_linear.pdf", height=5)
plot(owid.fit, inc="interval"
   , log=FALSE # linear scale
   ##, ylim=c(10^2,10^6)
   , annotate=FALSE
   , ylab="Daily reported cases"
   , line_style=list(col=0)
   , polygon_style=list(col=0)
     ##, main="Worldwide COVID-19 confirmed cases"
   , main=paste0("Worldwide COVID-19 confirmed cases up to ", get_date(owid.fit))
)
if (!interactive()) dev.off()

world_fits <- function(log=TRUE,
                       text_dbl1=list(y=200),
                       text_dbl4=list(x=245,y=10^5),
                       style.pt=list(points_main=list(bg="darkred",col="black")),
                       ...) {
## first wave:
owid.init1 <- update(owid.init, peak=38)
owid.fit1 <- egf(owid.init1)
## interval incidence
    ##plot(owid.init1, inc="interval")
    plot(owid.fit1, inc="interval"
       , log=log
       , ylim=if(log) c(10^2,10^6) else c(0,7e5)
       , annotate=FALSE
         ##, main="Worldwide COVID-19 confirmed cases"
       , main=paste0("Worldwide COVID-19 confirmed cases up to ", get_date(owid.fit))
       , ylab="Daily reported cases"
       , style=list(text_dbl=text_dbl1,
                    points_main=style.pt[["points_main"]])
       , ... )
    ##show_doubling_time(owid.fit1, pos=4, h=0.05, xshift=-15, cex=0.8)

## second wave:
owid.init2 <- update(owid.init, first=55, peak=90)
owid.fit2 <- egf(owid.init2)
## interval incidence
plot(owid.fit2, inc="interval", add=TRUE, style=style.pt)
    
##show_doubling_time(owid.fit2, pos=4, h=0.2, xshift=0, cex=0.8)

## third wave:
owid.init3 <- update(owid.init, first=130, peak=200)
owid.fit3 <- egf(owid.init3)
## interval incidence
plot(owid.fit3, inc="interval", add=TRUE, style=style.pt)
##show_doubling_time(owid.fit3, pos=4, h=0.1, xshift=-30, cex=0.8)

## fourth wave:
##owid.init4 <- update(owid.init, first=255)
owid.init4 <- update(owid.init, first=275)
owid.fit4 <- egf(owid.init4)
## interval incidence
    plot(owid.fit4, inc="interval", add=TRUE
       , style=list(text_dbl=text_dbl4,
                    points_main=style.pt[["points_main"]])
         )
##show_doubling_time(owid.fit4, pos=4, h=-0.1, xshift=-20, cex=0.8)

## show_date(owid.fit4)
}

## ALL STATUS PLOTS
if (!interactive()) {
    pdf("covid_status.pdf", height=5)
    world_fits(log=FALSE)
    world_fits(add=FALSE)
    ontario_fits(log=FALSE)
    ontario_fits()
    dev.off()
}
