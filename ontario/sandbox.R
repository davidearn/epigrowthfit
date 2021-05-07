library("epigrowthfit")
load("ontario.RData")
load("endpoints.RData")
options(contrasts = c("contr.sum", "contr.poly"))

object <- egf(
  formula_ts  = cases_new ~ Date | region,
  data_ts     = ontario,
  formula_par = ~I(region:window),
  data_par    = endpoints,
  endpoints   = endpoints,
  curve       = "logistic",
  distr       = "nbinom",
  debug       = TRUE
)
## FIXME:
## Below construction of `init` assumes that
## match(levels(endpoints$region), levels(ontario$region), 0L)
## is increasing
f <- function(x) c(m <- mean(x), x[-length(x)] - m)
init <- unlist(lapply(object$Y_init, f), use.names = FALSE)
object <- update(object, debug = FALSE, init = init)

pdf("ontario_incidence.pdf", width = 8, height = 4)
plot(object,
  type = "interval",
  show_tdoubling = TRUE,
  xlim = c("2020-03-01", "2021-05-01"),
  log = TRUE,
  order = order(population, decreasing = TRUE)
)
dev.off()
