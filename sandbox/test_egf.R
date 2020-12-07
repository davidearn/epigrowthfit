library("TMB")
setwd("src")
file.remove(paste0("epigrowthfit.", c("o", "so")))
compile("epigrowthfit.cpp")
dyn.load(dynlib("epigrowthfit"))
setwd("..")


load("data/canadacovid.RData")
ontario <- subset(canadacovid, province == "ON")
ontario$wave <- 0L
ontario$wave[27:77] <- 1L
ontario$wave[211:236] <- 2L
ontario$wave <- factor(ontario$wave, exclude = 0L)

for (filename in list.files("R", full.names = TRUE)) {
  source(filename)
}

object <- egf(new_confirmed ~ date,
  data = ontario,
  index = ontario$wave,
  #fixed = ~1,
  fixed = ~wave,
  #random = ~(1 | wave),
  random = NULL,
  spX = FALSE,
  spZ = FALSE,
  curve = "logistic",
  distr = "nbinom",
  include_baseline = FALSE,
  na_action = "fail",
  method = "nlminb",
  dfmt = "%Y-%m-%d"
)
coef(object)
coef(object, wave = 1)
vcov(object)
vcov(object, standardize = FALSE)
confint(object, parm = "doubling_time", level = 0.95, trace = FALSE)
predict(object, time = 0:20, se = TRUE, wave = 1)


