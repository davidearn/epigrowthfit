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

source("R/utils.R")
source("R/egf_checks.R")
source("R/egf_utils.R")
source("R/egf.R")
source("R/coef.R")
source("R/vcov.R")
source("R/predict.R")
source("R/profile.R")

object <- egf(new_confirmed ~ date,
  data = ontario,
  index = ontario$wave,
  fixed = ~1,
  #fixed = ~wave,
  random = ~(1 | wave),
  #random = NULL,
  sparse_X = FALSE,
  sparse_Z = FALSE,
  curve = "logistic",
  distr = "nbinom",
  excess = FALSE,
  na_action = "exclude",
  method = "nlminb",
  date_format = "%Y-%m-%d"
)

