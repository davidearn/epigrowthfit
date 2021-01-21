library("epigrowthfit")

ontario <- subset(canadacovid, province == "ON")
ontario$wave <- factor(rep(1:3, diff(c(0L, 100L, 250L, 346L))))
index <- rep(NA, nrow(ontario))
index[27:77] <- 1
index[211:245] <- 2
index[251:346] <- 3
index_ontario <- factor(index, exclude = NA)

quebec <- subset(canadacovid, province == "QC")
quebec$wave <- factor(rep(1:5, diff(c(0L, 100L, 150L, 220L, 255L, 313L))))
index <- rep(0, nrow(quebec))
index[1:35] <- 4
index[119:138] <- 5
index[165:215] <- 6
index[233:253] <- 7
index[260:310] <- 8
index_quebec <- factor(index, exclude = 0)

data <- droplevels(rbind(ontario, quebec))
index <- unlist(list(index_ontario, index_quebec))

object <- egf(new_confirmed ~ date,
  data = data,
  index = index,
  #fixed = ~1,
  fixed = ~province:wave,
  #random = ~(1 | wave),
  random = NULL,
  curve = "logistic",
  distr = "nbinom",
)

plot(object,
  group_by = ~province,
  control = list(text_dbl = list(cex = 0.6)),
  main = "Fitted logistic model (%province)",
  subset = list(province = "ON")
)

p <- profile(object, parm = "r")
