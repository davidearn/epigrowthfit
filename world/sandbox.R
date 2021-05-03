library("epigrowthfit")
load("data/world.RData")
load("data/endpoints.RData")
options(contrasts = c("contr.sum", "contr.poly"))
## FIXME:
## Below construction of `init` assumes that
## match(levels(endpoints$country_iso_alpha3), levels(world$country_iso_alpha3), 0L)
## is increasing
f <- function(x) {
  m <- mean(x)
  c(m, x[-length(x)] - m)
}
#endpoints <- droplevels(endpoints[complete.cases(endpoints), , drop = FALSE])

object <- egf(
  formula_ts  = cases_new ~ Date | country_iso_alpha3,
  data_ts     = world,
  formula_par = ~I(country_iso_alpha3:window),
  data_par    = endpoints,
  endpoints   = endpoints,
  curve       = "logistic",
  distr       = "nbinom",
  na_action   = c("exclude", "fail"),
  debug       = TRUE
)
#init <- rep_len(0, length(object_bigmod$init))
#i <- grepl("Intercept|country_iso_alpha3", object_bigmod$tmb_args$data$X_info$term)
#init[i] <- unlist(lapply(object_bigmod$Y_init, f), use.names = FALSE)
init <- unlist(lapply(object$Y_init, f), use.names = FALSE)
elapsed <- Sys.time()
object <- update(object, debug = FALSE, init = init)
elapsed <- Sys.time() - elapsed
save(object, elapsed, file = "object.RData")

pdf("world_incidence.pdf", width = 8, height = 4, onefile = TRUE)
plot(object_bigmod,
  type = "interval",
  show_tdoubling = TRUE,
  xlim = c("2020-01-01", "2021-05-01"),
  log = TRUE,
  sub = country
)
dev.off()


# pdf("world_heat_map_by_latitude.pdf", width = 6, height = 4, onefile = TRUE)
# plot(object,
#   type = "rt2",
#   per_plot = 15L,
#   xlim = c("2020-01-01", "2021-05-01"),
#   log = TRUE,
#   plab = country,
#   main = "Per capita growth rate, by country"
# )
# dev.off()


