library("epigrowthfit")
load("on_by_region.RData")
load("endpoints.RData")
options(contrasts = c("contr.sum", "contr.poly"))

window <- make_window(
  time = on_by_region$date,
  ts = on_by_region$region,
  endpoints = endpoints
)
wave <- make_wave(
  window = window,
  ts = on_by_region$region
)
on_by_region$`region:wave` <- (on_by_region$region):(wave)

object <- egf(
  formula_ts  = cases ~ date | region,
  formula_par = ~`region:wave`,
  data        = on_by_region,
  window      = window,
  curve       = "logistic",
  distr       = "nbinom",
  debug       = TRUE
)

f <- function(x) {m <- mean(x); c(m, x[-length(x)] - m)}
init <- unlist(lapply(object$Y_init, f), use.names = FALSE)
object <- update(object, debug = FALSE, init = init)

plot(object,
  type = "interval",
  show_tdoubling = TRUE,
  xlim = c("2020-03-01", "2021-03-15")
)
