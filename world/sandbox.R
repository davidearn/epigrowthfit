library("epigrowthfit")
load("world.RData")
load("endpoints.RData")
options(contrasts = c("contr.sum", "contr.poly"))

## For now, delete spurious zeros that would likely cause issues
## when estimating a model without weekday effects
omit_zeros <- function(d, origin) {
  zero <- d$date > as.Date(origin) & !is.na(d$cases) & d$cases == 0
  d[!zero, , drop = FALSE]
}
world_split <- split(world, world$country_iso3)
world_split$ESP <- omit_zeros(world_split$ESP, "2020-05-01")
world_split$BEL <- omit_zeros(world_split$BEL, "2020-05-01")
world_split$CZE <- omit_zeros(world_split$CZE, "2020-05-01")
world_split$SWE <- omit_zeros(world_split$SWE, "2020-05-01")
world_split$BLR <- omit_zeros(world_split$BLR, "2020-05-01")
world_split$CHE <- omit_zeros(world_split$CHE, "2020-05-01")
world_split$BGR <- omit_zeros(world_split$BGR, "2020-05-01")
world_split$SRB <- omit_zeros(world_split$SRB, "2020-04-01")
world_split$DNK <- omit_zeros(world_split$DNK, "2020-05-01")
world_split$FIN <- omit_zeros(world_split$FIN, "2020-11-01")
world_split$NOR <- omit_zeros(world_split$NOR, "2020-11-01")
world_split$BIH <- omit_zeros(world_split$BIH, "2020-06-01")
world_split$LUX <- omit_zeros(world_split$LUX, "2020-08-01")
world_split$AND <- omit_zeros(world_split$AND, "2020-08-01")
world_split$SMR <- omit_zeros(world_split$SMR, "2020-10-01")
world_split$MEX <- omit_zeros(world_split$MEX, "2021-01-01")
world_split$GTM <- omit_zeros(world_split$GTM, "2020-11-01")
world_split$HTI <- omit_zeros(world_split$HTI, "2020-06-01")
world_split$DOM <- omit_zeros(world_split$DOM, "2020-09-01")
world_split$HND <- omit_zeros(world_split$HND, "2020-05-01")
world_split$NIC <- omit_zeros(world_split$NIC, "2020-03-01")
world_split$SLV <- omit_zeros(world_split$SLV, "2020-09-01")
world_split$CRI <- omit_zeros(world_split$CRI, "2020-08-01")
world_split$BLZ <- omit_zeros(world_split$BLZ, "2020-09-01")
world_split$BHS <- omit_zeros(world_split$BHS, "2020-09-01")
world_split$BRA <- omit_zeros(world_split$BRA, "2020-09-01")
world_split$ARG <- omit_zeros(world_split$ARG, "2020-12-01")
world_split$PER <- omit_zeros(world_split$PER, "2020-06-01")
world_split$VEN <- omit_zeros(world_split$VEN, "2020-10-01")
world_split$ECU <- omit_zeros(world_split$ECU, "2020-04-01")
world_split$PNG <- omit_zeros(world_split$PNG, "2020-03-01")
world_split$IND <- omit_zeros(world_split$IND, "2020-12-01")
world_split$PAK <- omit_zeros(world_split$PAK, "2020-05-01")
world_split$TUR <- omit_zeros(world_split$TUR, "2020-12-01")
world_split$THA <- omit_zeros(world_split$THA, "2021-01-01")
world_split$UZB <- omit_zeros(world_split$UZB, "2020-11-01")
world <- droplevels(do.call(rbind, world_split))
row.names(world) <- NULL

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
f <- function(x) {m <- mean(x); c(m, x[-length(x)] - m)}
init <- unlist(lapply(object$Y_init, f), use.names = FALSE)
object <- update(object, debug = FALSE, init = init)
save(object, file = "object.RData")
load("object.RData")

pdf("world_incidence.pdf", width = 8, height = 4, onefile = TRUE)
plot(object,
  type = "interval",
  show_tdoubling = TRUE,
  xlim = c("2020-01-01", "2021-04-01"),
  log = TRUE,
  order = order(population, decreasing = TRUE),
  sub = country_name
)
dev.off()

pdf("world_heat_map.pdf", width = 6, height = 4, onefile = TRUE)
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
