library("epigrowthfit")
load("world.RData")

europe <- droplevels(subset(world, continent == "Europe"))
europe <- europe[names(europe) != "continent"]
row.names(europe) <- NULL

## Choose fitting windows by eye
europe_split <- split(europe, europe$country_name)
dline <- function(date0, date, ...) {
  i <- which.min(abs(date - as.Date(date0)))
  abline(v = date[i], ...)
}

s <- "Italy"
d <- europe_split[[s]]
i <- is.finite(d$cases)
ss <- smooth.spline(
  x = d$date[i] - d$date[1L],
  y = log10(1 + d$cases[i]),
  spar = 0.5
)
plot(x = d$date, y = log10(1 + d$cases), pch = 16, col = "#66666680", main = s, xaxt = "n")
lines(x = d$date[1L] + ss$x, y = ss$y, lwd = 2, col = "red")
at <- seq(as.Date("2020-01-01"), as.Date("2021-02-01"), by = "1 month")
axis(side = 1, at = at, labels = rep_len(month.abb, 14L))
abline(v = at, lty = 3, lwd = 1, col = "grey70")
dline("2020-12-15", date = d$date, col = "blue")
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
    c("2020-12-05", "2020-01-10")
  )
)
