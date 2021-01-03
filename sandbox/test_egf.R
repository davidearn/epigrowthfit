load("data/canadacovid.RData")
ontario <- subset(canadacovid, province == "ON")
ontario$wave <- 0L
ontario$wave[27:77] <- 1L
ontario$wave[211:236] <- 2L
ontario$wave <- factor(ontario$wave)

# source("R/predict.R")
# source("R/profile.R")
# source("R/confint.R")
# source("R/boot.R")

object <- egf(new_confirmed ~ date,
  data = ontario,
  index = factor(ontario$wave, exclude = 0L),
  #fixed = ~1,
  fixed = ~wave,
  #random = ~(1 | wave),
  random = NULL,
  curve = "logistic",
  distr = "nbinom",
)

