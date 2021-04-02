library("epigrowthfit")
load("on_by_region.RData")
load("endpoints.RData")
options(contrasts = c("contr.sum", "contr.poly"))

endpoints$`region:wave` <- endpoints$region:endpoints$wave

object <- egf(
  formula_ts  = cases ~ date | region,
  data_ts     = on_by_region,
  formula_par = ~`region:wave`,
  data_par    = endpoints,
  endpoints   = endpoints,
  curve       = "logistic",
  distr       = "nbinom",
  debug       = TRUE
)

f <- function(x) {m <- mean(x); c(m, x[-length(x)] - m)}
init <- unlist(lapply(object$Y_init, f), use.names = FALSE)
object <- update(object, debug = FALSE, init = init, append = c(region, wave))

plot(object,
  type = "interval",
  show_tdoubling=TRUE,
  xlim = c("2020-03-01", "2021-03-15"),
  log=TRUE
)


plot(object,
  type = "rt2",
  per_plot = 12L,
  xlim = c("2020-03-01", "2021-03-15"),
  log=FALSE
)
