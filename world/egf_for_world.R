library("epigrowthfit")
load("world.RData")
europe <- droplevels(subset(world, continent == "Europe"))
europe <- europe[names(europe) != "continent"]

## Omit unexpected zeros
countries <- c(
  "Spain",
  "Belgium",
  "Czech Republic",
  "Sweden",
  "Belarus",
  "Switzerland",
  "Bulgaria",
  "Serbia",
  "Finland",
  "Norway"
)
starts <- as.Date(
  rep.int(c("2020-05-01", "2020-11-01"), c(8L, 2L))
)
is_zero_cases <- (!is.na(europe$cases) & europe$cases == 0L)
f <- function(country, start) europe$country == country & europe$date >= start & is_zero_cases
europe <- europe[!Reduce(`|`, Map(f, countries, starts)), , drop = FALSE]
row.names(europe) <- NULL

## Choose fitting windows by eye
# europe_split <- split(europe, europe$country)
# dline <- function(date0, date, ...) {
#   i <- vapply(date0, function(x) which.min(abs(date - as.Date(x))), 0L)
#   abline(v = date[i], ...)
# }
# do_plot <- function(s, spar) {
#   d <- europe_split[[s]]
#   i <- is.finite(d$cases)
#   ss <- smooth.spline(
#     x = d$date[i] - d$date[1L],
#     y = log10(1 + d$cases[i]),
#     spar = spar
#   )
#   plot(x = d$date, y = log10(1 + d$cases), pch = 16, col = "#66666680", main = s, xaxt = "n", las = 1, mgp = c(2.5, 0.7, 0), xlab = "", ylab = expression(log[10](1 + cases)))
#   lines(x = d$date[1L] + ss$x, y = ss$y, lwd = 2, col = "red")
#   at <- seq(as.Date("2020-01-01"), as.Date("2021-02-01"), by = "1 month")
#   axis(side = 1, at = at, labels = rep_len(month.abb, 14L))
#   abline(v = at, lty = 3, lwd = 1, col = "grey70")
# }
#
# country <- "Norway"
# spar <- 0.5
# do_plot(country, spar)
# dline(c("2020-12-12", "2021-01-10"), date = europe_split[[country]]$date, col = "blue")

w <- list(
  Russia = list(
    c("2020-03-05", "2020-05-10"),
    c("2020-08-25", "2020-09-20"),
    c("2020-09-22", "2020-10-23"),
    c("2020-10-28", "2020-12-15")
  ),
  Germany = list(
    c("2020-02-25", "2020-04-05"),
    c("2020-07-10", "2020-08-23"),
    c("2020-09-01", "2020-09-30"),
    c("2020-10-01", "2020-11-10"),
    c("2020-11-25", "2020-12-25")
  ),
  France = list(
    c("2020-02-25", "2020-04-05"),
    c("2020-07-01", "2020-09-27"),
    c("2020-10-01", "2020-11-10"),
    c("2020-11-28", "2021-01-15")
  ),
  `United Kingdom` = list(
    c("2020-02-20", "2020-04-10"),
    c("2020-07-05", "2020-08-28"),
    c("2020-08-31", "2020-11-01"),
    c("2020-12-05", "2021-01-10")
  ),
  Italy = list(
    c("2020-02-23", "2020-03-25"),
    c("2020-08-05", "2020-09-07"),
    c("2020-09-28", "2020-11-15")
  ),
  Spain = list(
    c("2020-02-25", "2020-04-01"),
    c("2020-06-26", "2020-10-01"),
    c("2020-10-10", "2020-11-07"),
    c("2020-12-27", "2021-01-25")
  ),
  Ukraine = list(
    c("2020-03-15", "2020-04-28"),
    c("2020-05-24", "2020-06-29"),
    c("2020-07-15", "2020-08-24"),
    c("2020-08-25", "2020-12-01")
  ),
  Poland = list(
    c("2020-03-01", "2020-04-10"),
    c("2020-07-13", "2020-08-13"),
    c("2020-09-06", "2020-11-10")
  ),
  Romania = list(
    c("2020-03-05", "2020-04-13"),
    c("2020-05-28", "2020-06-30"),
    c("2020-07-01", "2020-08-01"),
    c("2020-09-15", "2020-11-13")
  ),
  Netherlands = list(
    c("2020-02-25", "2020-04-01"),
    c("2020-07-06", "2020-08-15"),
    c("2020-08-28", "2020-11-01"),
    c("2020-11-30", "2020-12-23")
  ),
  Belgium = list(
    c("2020-03-01", "2020-04-06"),
    c("2020-07-01", "2020-08-12"),
    c("2020-09-01", "2020-09-28"),
    c("2020-09-29", "2020-11-01")
  ),
  Greece = list(
    c("2020-02-25", "2020-04-03"),
    c("2020-06-25", "2020-07-15"),
    c("2020-07-20", "2020-08-18"),
    c("2020-09-01", "2020-09-30"),
    c("2020-10-10", "2020-11-15"),
    c("2021-01-18", "2021-02-04")
  ),
  `Czech Republic` = list(
    c("2020-03-01", "2020-04-03"),
    c("2020-06-10", "2020-07-01"),
    c("2020-07-10", "2020-07-29"),
    c("2020-08-20", "2020-09-25"),
    c("2020-09-28", "2020-10-30"),
    c("2020-12-01", "2021-01-10")
  ),
  Sweden = list(
    c("2020-02-25", "2020-04-12"),
    c("2020-05-25", "2020-06-20"),
    c("2020-07-20", "2020-08-12"),
    c("2020-08-31", "2020-10-15"),
    c("2020-10-16", "2020-11-15"),
    c("2020-11-16", "2021-01-08")
  ),
  Portugal = list(
    c("2020-03-01", "2020-04-03"),
    c("2020-05-05", "2020-06-12"),
    c("2020-08-20", "2020-09-20"),
    c("2020-10-01", "2020-11-19"),
    c("2020-12-24", "2021-01-30")
  ),
  Hungary = list(
    c("2020-03-01", "2020-04-15"),
    c("2020-06-15", "2020-08-18"),
    c("2020-08-20", "2020-09-22"),
    c("2020-10-12", "2020-11-18")
  ),
  Belarus = list(
    c("2020-04-01", "2020-05-01"),
    c("2020-08-13", "2020-08-28"),
    c("2020-09-15", "2020-11-01"),
    c("2020-11-08", "2020-12-13")
  ),
  Austria = list(
    c("2020-02-20", "2020-03-28"),
    c("2020-06-15", "2020-07-20"),
    c("2020-08-01", "2020-08-23"),
    c("2020-09-01", "2020-09-23"),
    c("2020-10-12", "2020-11-14")
  ),
  Switzerland = list(
    c("2020-02-20", "2020-03-28"),
    c("2020-06-15", "2020-07-10"),
    c("2020-07-12", "2020-09-18"),
    c("2020-09-25", "2020-11-06")
  ),
  Bulgaria = list(
    c("2020-03-01", "2020-04-01"),
    c("2020-04-10", "2020-04-28"),
    c("2020-05-25", "2020-06-15"),
    c("2020-06-16", "2020-07-23"),
    c("2020-09-20", "2020-11-15"),
    c("2021-01-17", "2021-02-05")
  ),
  Serbia = list(
    c("2020-03-10", "2020-04-17"),
    c("2020-06-15", "2020-07-10"),
    c("2020-10-15", "2020-12-01")
  ),
  Denmark = list(
    c("2020-07-07", "2020-08-15"),
    c("2020-08-25", "2020-09-27"),
    c("2020-10-17", "2020-11-07"),
    c("2020-11-27", "2020-12-20")
  ),
  Finland = list(
    c("2020-03-01", "2020-04-10"),
    c("2020-07-10", "2020-10-18"),
    c("2021-01-12", "2021-02-05")
  ),
  Slovakia = list(
    c("2020-03-01", "2020-04-24"),
    c("2020-06-01", "2020-07-12"),
    c("2020-07-15", "2020-09-12"),
    c("2020-09-15", "2020-11-01"),
    c("2020-12-01", "2020-12-20")
  ),
  Norway = list(
    c("2020-02-25", "2020-03-30"),
    c("2020-07-16", "2020-08-13"),
    c("2020-08-25", "2020-09-11"),
    c("2020-10-18", "2020-11-15"),
    c("2020-12-12", "2021-01-10")
  )
)

ff <- make_index(
  date = europe$date,
  ts = europe$country,
  subset = w
)
europe$wave <- ff$wave
index <- ff$index

object <- egf(
  formula  = cases ~ date,
  fixed    = ~country:wave,
  random   = NULL,
  group_by = ~country,
  data     = europe,
  index    = index,
  curve    = "logistic",
  distr    = "nbinom"
)

ci <- confint(object, parm = "tdoubling", link = FALSE)
endpoints <- epigrowthfit:::make_endpoints(object)
endpoints[] <- lapply(endpoints, `+`, attr(endpoints, "refdate"))
europe_td <- cbind(ci, `names<-`(endpoints, c("date1", "date2")))
save(europe_td, file = "europe_td.RData")

pdf("europe_td1.pdf", width = 6, height = 4, onefile = TRUE)
plot(ci, type = "1", group_by = ~country, per_plot = 12)
dev.off()

pdf("europe_td2.pdf", width = 6, height = 4, onefile = TRUE)
plot(ci, type = "1", group_by = ~wave, sort = "increasing", per_plot = 12)
dev.off()

subset <- list(country = names(w))
xlim <- c("2020-02-05", "2021-02-05")

pdf("europe_incidence.pdf", width = 8, height = 4, onefile = TRUE)
plot(object,
  type = "interval",
  bands = TRUE,
  subset = subset,
  tol = Inf,
  xlim = xlim,
  ylab = "reported incidence",
  main = "Daily reported COVID-19 incidence\n%country"
)
dev.off()

pdf("on_rt1.pdf", width = 8, height = 4, onefile = TRUE)
plot(object,
  type = "rt1",
  bands = TRUE,
  subset = subset,
  xlim = xlim,
  ylim = c(0, log(2) / 2.5),
  main = "Instantaneous exponential growth rate\n%country",
  control = list(points = NULL, abline = NULL)
)
dev.off()

pdf("europe_rt2.pdf", width = 6, height = 4, onefile = TRUE)
plot(object,
  type = "rt2",
  per_plot = 25,
  bw_panels = 0,
  subset = subset,
  control = list(colorRamp = list(bias = 2))
)
dev.off()



