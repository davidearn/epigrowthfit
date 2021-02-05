library("epigrowthfit")
load("../on_by_region.RData")

## Ordered by population size
census_divisions <- c("TORONTO", "PEEL", "YORK", "OTTAWA",
                      "DURHAM", "HALTON", "HAMILTON", "WATERLOO",
                      "SIMCOE", "MIDDLESEX", "NIAGARA", "ESSEX")
w <- list(
  TORONTO = list(
    c("2020-03-05", "2020-04-23"),
    c("2020-08-07", "2020-09-05"),
    c("2020-09-06", "2020-10-05"),
    c("2020-10-11", "2021-01-13")
  ),
  PEEL = list(
    c("2020-03-14", "2020-04-23"),
    c("2020-08-02", "2020-10-10"),
    c("2020-10-16", "2020-11-19"),
    c("2020-12-25", "2021-01-08")
  ),
  YORK = list(
    c("2020-03-09", "2020-04-18"),
    c("2020-08-12", "2020-10-26"),
    c("2020-11-01", "2020-11-20"),
    c("2020-12-11", "2021-01-09")
  ),
  OTTAWA = list(
    c("2020-03-09", "2020-04-25"),
    c("2020-07-01", "2020-07-31"),
    c("2020-08-31", "2020-10-09"),
    c("2020-12-19", "2021-01-17")
  ),
  DURHAM = list(
    c("2020-03-17", "2020-04-25"),
    c("2020-08-16", "2020-10-25"),
    c("2020-10-31", "2021-01-13")
  ),
  HALTON = list(
    c("2020-03-16", "2020-04-10"),
    c("2020-08-30", "2020-10-19"),
    c("2020-10-25", "2020-11-13"),
    c("2020-11-24", "2021-01-07")
  ),
  HAMILTON = list(
    c("2020-03-10", "2020-04-15"),
    c("2020-09-09", "2020-10-23"),
    c("2020-10-29", "2021-01-06")
  ),
  WATERLOO = list(
    c("2020-03-03", "2020-04-25"),
    c("2020-09-04", "2020-09-28"),
    c("2020-10-24", "2020-11-27"),
    c("2020-12-03", "2021-01-11")
  ),
  SIMCOE = list(
    c("2020-03-11", "2020-04-17"),
    c("2020-08-20", "2020-10-12"),
    c("2020-10-18", "2021-01-15")
  ),
  MIDDLESEX = list(
    c("2020-03-16", "2020-04-14"),
    c("2020-09-07", "2020-10-16"),
    c("2020-11-01", "2021-01-14")
  ),
  NIAGARA = list(
    c("2020-03-12", "2020-04-16"),
    c("2020-09-06", "2020-11-15"),
    c("2020-11-26", "2021-01-14")
  ),
  ESSEX = list(
    c("2020-03-16", "2020-04-12"),
    c("2020-10-20", "2021-01-02")
  )
)

## Choose fitting windows by eye
# on_by_region_split <- split(on_by_region, factor(on_by_region$region, levels = census_divisions))
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
on_td <- cbind(ci, `names<-`(endpoints, c("date1", "date2")))
save(on_td, file = "on_td.RData")

pdf("on_td1.pdf", width = 6, height = 4, onefile = TRUE)
plot(ci, type = "1", group_by = ~region, per_plot = 12)
dev.off()

pdf("on_td2.pdf", width = 6, height = 4, onefile = TRUE)
plot(ci, type = "1", group_by = ~wave, sort = "increasing", per_plot = 12)
dev.off()

subset <- list(region = census_divisions)
xlim <- c("2020-03-01", "2021-02-02")

pdf("on_incidence.pdf", width = 8, height = 4, onefile = TRUE)
plot(object,
  type = "interval",
  bands = TRUE,
  subset = subset,
  tol = Inf,
  xlim = xlim,
  ylab = "reported incidence",
  main = "Daily reported COVID-19 incidence\n%region"
)
dev.off()

pdf("on_rt1.pdf", width = 8, height = 4, onefile = TRUE)
plot(object,
  type = "rt1",
  bands = TRUE,
  subset = subset,
  xlim = xlim,
  ylim = c(0, log(2) / 2.5),
  main = "Instantaneous exponential growth rate\n%region",
  control = list(points = NULL, abline = NULL)
)
dev.off()

pdf("on_rt2.pdf", width = 6, height = 4, onefile = TRUE)
plot(object,
  type = "rt2",
  per_plot = 12,
  bw_panels = 0,
  subset = subset,
  control = list(colorRamp = list(bias = 2))
)
dev.off()
