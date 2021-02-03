library("epigrowthfit")
load("on_by_region.RData")

## Ordered by population size
census_divisions <- c("TORONTO", "PEEL", "YORK", "OTTAWA",
                      "DURHAM", "HALTON", "HAMILTON", "WATERLOO",
                      "SIMCOE", "MIDDLESEX", "NIAGARA", "ESSEX")
on_by_region_split <- split(on_by_region, factor(on_by_region$region, levels = census_divisions))
w <- list(
  TORONTO = list(11:60, 166:195, 196:225, 231:325),
  PEEL = list(6:45, 146:215, 221:255, 291:305),
  YORK = list(6:45, 161:235, 241:260, 281:310),
  OTTAWA = list(1:45, 111:140, 171:210, 281:310),
  DURHAM = list(6:45, 156:225, 231:305),
  HALTON = list(6:30, 161:210, 216:235, 246:290),
  HAMILTON = list(1:35, 171:215, 221:290),
  WATERLOO = list(1:45, 161:185, 211:245, 251:290),
  SIMCOE = list(1:35, 141:190, 196:285),
  MIDDLESEX = list(3:30, 151:190, 206:280),
  NIAGARA = list(1:30, 151:220, 231:280),
  ESSEX = list(1:25, 216:290)
)

## Choose fitting windows by eye
# s <- "TORONTO"
# d <- on_by_region_split[[s]]
# ss <- smooth.spline(
#   x = d$date - d$date[1],
#   y = log10(1 + d$cases),
#   spar = 0.7
# )
# plot(x = d$date - d$date[1], y = log10(1 + d$cases), main = s)
# lines(ss, lwd = 2, col = "red")
# i <- 1
# abline(v = d$date[i] - d$date[1])

ff <- make_index(
  date = on_by_region$date,
  ts = on_by_region$region,
  subset = w
)
on_by_region$wave <- ff$wave
index <- ff$index

object <- egf(
  formula  = cases ~ date,
  fixed    = ~region:wave,
  random   = NULL,
  group_by = ~region,
  data     = on_by_region,
  index    = index,
  curve    = "logistic",
  distr    = "nbinom"
)

ci <- confint(object, parm = "tdoubling", link = FALSE)
endpoints <- epigrowthfit:::make_endpoints(object)
endpoints[] <- lapply(endpoints, `+`, attr(endpoints, "refdate"))
ontario_doubling_times <- cbind(ci, `names<-`(endpoints, c("date1", "date2")))
save(ontario_doubling_times, file = "ontario_doubling_times.RData")

pdf("ontario_doubling1.pdf", width = 6, height = 4, onefile = TRUE)
plot(ci, type = "1", group_by = ~region, per_plot = 12)
dev.off()

pdf("ontario_doubling2.pdf", width = 6, height = 6, onefile = TRUE)
plot(ci, type = "2", per_plot = 6)
dev.off()

subset <- list(region = census_divisions)
subset <- list(region = "TORONTO")

pdf("ontario_incidence.pdf", width = 8, height = 4, onefile = TRUE)
plot(object,
  type = "interval",
  bands = TRUE,
  subset = subset,
  tol = Inf
)
dev.off()

pdf("ontario_rt1.pdf", width = 8, height = 4, onefile = TRUE)
plot(object,
  type = "rt1",
  bands =  TRUE,
  subset = subset,
  ylim = c(0, 0.2),
  control = list(points = NULL)
)
dev.off()

pdf("ontario_rt2.pdf", width = 6, height = 4, onefile = TRUE)
plot(object,
  type = "rt2",
  per_plot = 12,
  bw_plot = 0,
  subset = subset,
  control = list(heat = list(bias = 2))
)
