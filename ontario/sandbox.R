library("epigrowthfit")
load("ontario.RData")
load("endpoints.RData")
options(contrasts = c("contr.sum", "contr.poly"))

object <- egf(
  formula_ts  = cases ~ date | region,
  data_ts     = ontario,
  formula_par = ~I(region:wave),
  data_par    = endpoints,
  endpoints   = endpoints,
  curve       = "logistic",
  distr       = "nbinom",
  append      = -c(start, end),
  debug       = TRUE
)
## FIXME:
## Below construction of `init` assumes that
## match(levels(endpoints$region), levels(ontario$region), 0L)
## is increasing
f <- function(x) c(m <- mean(x), x[-length(x)] - m)
init <- unlist(lapply(object$Y_init, f), use.names = FALSE)
object <- update(object, debug = FALSE, init = init)

plot(object,
  type = "interval",
  show_tdoubling = TRUE,
  xlim = c("2020-03-01", "2021-03-25"),
  log = TRUE,
  order = order(population, decreasing = TRUE)
)

plot(object,
  type = "rt2",
  per_plot = 12L,
  xlim = c("2020-03-01", "2021-03-25"),
  log = TRUE,
  order = order(population, decreasing = TRUE)
)
