library("epigrowthfit")
load("world.RData")
load("endpoints.RData")
options(contrasts = c("contr.sum", "contr.poly"))

object <- egf(
  formula_ts  = cases ~ date | country_iso3,
  data_ts     = world,
  formula_par = ~I(country_iso3:wave),
  data_par    = endpoints,
  endpoints   = endpoints,
  curve       = "logistic",
  distr       = "nbinom",
  append      = -c(start, end),
  na_action   = c("exclude", "fail"),
  debug       = TRUE
)
## FIXME:
## Below construction of `init` assumes that
## match(levels(endpoints$country_iso3), levels(world$country_iso3), 0L)
## is increasing
f <- function(x) {
  m <- mean(x)
  c(m, x[-length(x)] - m)
}
init <- unlist(lapply(object$Y_init, f), use.names = FALSE)
elapsed <- Sys.time()
object <- update(object, debug = FALSE, init = init)
elapsed <- Sys.time() - elapsed
save(object, elapsed, file = "object.RData")
# load("object.RData")

pdf("world_incidence.pdf", width = 8, height = 4, onefile = TRUE)
plot(object,
  type = "interval",
  show_tdoubling = TRUE,
  xlim = c("2020-01-01", "2021-04-10"),
  log = TRUE,
  order = order(latitude, decreasing = TRUE),
  sub = country_name
)
dev.off()

pdf("world_heat_map_by_latitude.pdf", width = 6, height = 4, onefile = TRUE)
plot(object,
  type = "rt2",
  per_plot = 15L,
  xlim = c("2020-01-01", "2021-04-01"),
  log = TRUE,
  order = order(latitude, decreasing = TRUE),
  plab = country_name,
  main = "Per capita growth rate, by country"
)
dev.off()
