library("epigrowthfit")

##if (!interactive()) pdf("covid_growth.pdf", height=5.5, width=8.5)

## A helper for fitting window selection ... example of usage below
do_plot <- function(data, spar = NULL, window = NULL, do.log=TRUE, ...) {
  op <- par(mar = c(2, 3, 1, 1), mgp = c(2, 0.7, 0), las = 1)
  on.exit(par(op))
  xat <- seq(as.Date("2020-01-01"), Sys.Date(), by = "1 month")
  if (do.log) {
      plot(log1p(cases) ~ date, data = data, xlab = "", xaxt = "n", ...)
  } else {
      plot(cases ~ date, data = data, xlab = "", xaxt = "n", ...)
  }
  abline(v = xat, lty = 3, col = "grey80")
  axis(side = 1, at = xat, labels = months(xat, TRUE), gap.axis = 0)
  if (!is.null(spar)) {
      dd <- na.omit(
          data.frame(x = data$date - data$date[1],
                     y = if (do.log) log1p(data$cases) else data$cases))
    ss <- smooth.spline(x = dd$x, y = dd$y, spar = spar)
    lines(x = data$date[1] + ss$x, y = ss$y, lwd = 2, col = "red")
  }
  if (!is.null(window)) {
    abline(v = as.Date(unlist(window)), lty = 3)
  }
  invisible(NULL)
}


## Ontario plot

canada <- read.csv("COVID19_Canada.csv")
ontario <- subset(canada, Province == "ON")
ontario <- data.frame(
  date = as.Date(ontario$Date),
  cases = c(NA, diff(ontario$confirmed_positive))
)
w <- list(
  c("2020-03-01", "2020-04-10"), # wave 1
  c("2020-09-01", "2020-10-02"), # wave 2
  c("2020-10-05", "2021-01-10")  # wave 3
)
do_plot(data = ontario, spar = 0.6, window = w, do.log=FALSE)
do_plot(data = ontario, spar = 0.6, window = w)

ff <- make_index(
  date = ontario$date,
  ts = rep(factor("ontario"), nrow(ontario)),
  subset = list(ontario = w)
)

egf_ontario <- egf(
  formula = cases ~ date,
  fixed = ~wave,
  data = cbind(ontario, wave = ff$wave),
  index = ff$index
)
if (!interactive()) pdf("covid_growth_ON_linear.pdf", height=5.5, width=8.5)
plot(egf_ontario, log=FALSE,
  bands = TRUE,
  ylab = "new cases",
  main = "Daily COVID-19 incidence\nOntario"
)
dev.off()
if (!interactive()) pdf("covid_growth_ON_log.pdf", height=5.5, width=8.5)
plot(egf_ontario,
  bands = TRUE,
  ylab = "new cases",
  main = "Daily COVID-19 incidence\nOntario"
)
dev.off()


## World plot ... once you have defined windows (`w`),
## commented code should run

world <- read.csv("owid-covid-data.csv")
tt <- tapply(world$new_cases, INDEX = world$date, FUN = sum, na.rm = TRUE)
world <- data.frame(
  date = as.Date(names(tt)),
  cases = unname(tt)
)

# w <- list(
#   c("2020-03-01", "2020-04-10"), # wave 1
#   c("2020-09-01", "2020-10-02"), # wave 2
#   c("2020-10-05", "2021-01-10")  # wave 3
# )
do_plot(data = world, spar = 0.6)
# do_plot(data = world, spar = 0.6, window = w)
#
# ff <- make_index(
#   date = world$date,
#   ts = rep(factor("world"), nrow(world)),
#   subset = list(world = w)
# )
#
# egf_world <- egf(
#   formula = cases ~ date,
#   fixed = ~wave,
#   data = cbind(world, wave = ff$wave),
#   index = ff$index
# )
# plot(egf_world,
#   bands = TRUE,
#   ylab = "new cases",
#   main = "Daily COVID-19 incidence\nWorld"
# )

##if (!interactive()) dev.off()
