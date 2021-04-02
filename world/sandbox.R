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
world_split <- split(world, world$country)
world_split$Spain                    <- omit_zeros(world_split$Spain,                    "2020-05-01")
world_split$Belgium                  <- omit_zeros(world_split$Belgium,                  "2020-05-01")
world_split$`Czech Republic`         <- omit_zeros(world_split$`Czech Republic`,         "2020-05-01")
world_split$Sweden                   <- omit_zeros(world_split$Sweden,                   "2020-05-01")
world_split$Belarus                  <- omit_zeros(world_split$Belarus,                  "2020-05-01")
world_split$Switzerland              <- omit_zeros(world_split$Switzerland,              "2020-05-01")
world_split$Bulgaria                 <- omit_zeros(world_split$Bulgaria,                 "2020-05-01")
world_split$Serbia                   <- omit_zeros(world_split$Serbia,                   "2020-04-01")
world_split$Denmark                  <- omit_zeros(world_split$Denmark,                  "2020-05-01")
world_split$Finland                  <- omit_zeros(world_split$Finland,                  "2020-11-01")
world_split$Norway                   <- omit_zeros(world_split$Norway,                   "2020-11-01")
world_split$`Bosnia and Herzegovina` <- omit_zeros(world_split$`Bosnia and Herzegovina`, "2020-06-01")
world_split$Luxembourg               <- omit_zeros(world_split$Luxembourg,               "2020-08-01")
world_split$Andorra                  <- omit_zeros(world_split$Andorra,                  "2020-08-01")
world_split$`San Marino`             <- omit_zeros(world_split$`San Marino`,             "2020-10-01")
world_split$Mexico                   <- omit_zeros(world_split$Mexico,                   "2021-01-01")
world_split$Guatemala                <- omit_zeros(world_split$Guatemala,                "2020-11-01")
world_split$Haiti                    <- omit_zeros(world_split$Haiti,                    "2020-06-01")
world_split$`Dominican Republic`     <- omit_zeros(world_split$`Dominican Republic`,     "2020-09-01")
world_split$Honduras                 <- omit_zeros(world_split$Honduras,                 "2020-05-01")
world_split$Nicaragua                <- omit_zeros(world_split$Nicaragua,                "2020-03-01")
world_split$`El Salvador`            <- omit_zeros(world_split$`El Salvador`,            "2020-09-01")
world_split$`Costa Rica`             <- omit_zeros(world_split$`Costa Rica`,             "2020-08-01")
world_split$Belize                   <- omit_zeros(world_split$Belize,                   "2020-09-01")
world_split$Bahamas                  <- omit_zeros(world_split$Bahamas,                  "2020-09-01")
world_split$Brazil                   <- omit_zeros(world_split$Brazil,                   "2020-09-01")
world_split$Argentina                <- omit_zeros(world_split$Argentina,                "2020-12-01")
world_split$Peru                     <- omit_zeros(world_split$Peru,                     "2020-06-01")
world_split$Venezuela                <- omit_zeros(world_split$Venezuela,                "2020-10-01")
world_split$Ecuador                  <- omit_zeros(world_split$Ecuador,                  "2020-04-01")
world_split$`Papua New Guinea`       <- omit_zeros(world_split$`Papua New Guinea`,       "2020-03-01")
world <- droplevels(do.call(rbind, world_split))
row.names(world) <- NULL

window <- make_window(
  time = world$date,
  ts = world$country,
  endpoints = endpoints
)
wave <- make_wave(
  window = window,
  ts = world$country
)
world$`country:wave` <- (world$country):(wave)


object <- egf(
  formula_ts  = cases ~ date | country,
  formula_par = ~`country:wave`,
  data        = world,
  window      = window,
  curve       = "logistic",
  distr       = "nbinom",
  debug       = TRUE
)

f <- function(x) {m <- mean(x); c(m, x[-length(x)] - m)}
init <- unlist(lapply(object$Y_init, f), use.names = FALSE)
object <- update(object, debug = FALSE, init = init)
save(object, file = "object.RData")
load("object.RData")

pdf("world_incidence.pdf", width = 8, height = 4, onefile = TRUE)
cache <- plot(object,
  type = "interval",
  show_tdoubling = TRUE,
  xlim = c("2020-02-25", "2021-03-25")
)
dev.off()


pdf("world_heat_map.pdf", width = 6, height = 4, onefile = TRUE)
plot(object,
  type = "rt2",
  per_plot = 12L,
  xlim = c("2020-03-01", "2021-03-15"),
  log = TRUE
)
dev.off()


confint(object)
