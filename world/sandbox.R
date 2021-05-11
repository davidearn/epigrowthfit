library("epigrowthfit")
load("data/world.RData")
load("data/endpoints.RData")
options(contrasts = c("contr.sum", "contr.poly"), warn = 1)
outfile <- file("sandbox.Rout", open = "wt")
sink(outfile, type = "output")
sink(outfile, type = "message")

mm <- grep("^mobility_", names(endpoints), value = TRUE)
cc <- complete.cases(endpoints[mm])
endpoints <- droplevels(endpoints[cc, , drop = FALSE])

## FIXME:
## Below construction of `init` assumes that
## match(levels(endpoints$country_iso_alpha3), levels(world$country_iso_alpha3), 0L)
## is increasing
f <- function(x) {
  m <- mean(x)
  c(m, x[-length(x)] - m)
}

if (file.exists("object.RData")) {
  load("object.RData")
} else {
  Y_init <- `dim<-`(numeric(0L), c(0L, 4L))
  for (d in split(endpoints, endpoints$country_iso_alpha3)) {
    object <- egf(
      formula_ts  = cases_new ~ Date,
      data_ts     = world,
      formula_par = if (nrow(d) > 1L) ~window else ~1,
      data_par    = d,
      endpoints   = d,
      curve       = "logistic",
      distr       = "nbinom",
      na_action   = c("exclude", "fail"),
      debug       = TRUE
    )
    init <- unlist(lapply(object$Y_init, f), use.names = FALSE)
    object <- update(object, debug = FALSE, init = init)
    Y_init <- rbind(Y_init, matrix(fitted(object)$estimate, ncol = 4L))
  }
  init <- unlist(lapply(as.data.frame(Y_init), f), use.names = FALSE)

  object <- egf(
    formula_ts  = cases_new ~ Date | country_iso_alpha3,
    data_ts     = world,
    formula_par = ~I(country_iso_alpha3:window),
    data_par    = endpoints,
    endpoints   = endpoints,
    curve       = "logistic",
    distr       = "nbinom",
    na_action   = c("exclude", "fail"),
    init        = init
  )
  save(object, file = "object.RData")

  pdf("world_incidence.pdf", width = 8, height = 4, onefile = TRUE)
  plot(object,
    type = "interval",
    show_tdoubling = TRUE,
    xlim = c(as.Date("2020-01-01"), Sys.Date()),
    log = TRUE,
    sub = country
  )
  dev.off()
}

if (!file.exists("object_m.RData")) {
  load("object_m.RData")
} else {
  Y_init <- matrix(fitted(object)$estimate, ncol = 4L)
  init <- c(f(Y_init[, 1L]), 0, 0, 0, 0, 0, 0, colMeans(Y_init[, -1L]), 1, 1, 1)
  object_m <- egf(
    formula_ts  = cases_new ~ Date | country_iso_alpha3,
    data_ts     = world,
    formula_par = list(
      log(r) ~
        I(country_iso_alpha3:window) +
        log(mobility_retail_and_recreation) +
        log(mobility_grocery_and_pharmacy) +
        log(mobility_parks) +
        log(mobility_transit_stations) +
        log(mobility_workplaces) +
        log(mobility_residential),
      log(tinfl)  ~ 1 | country_iso_alpha3:window,
      log(K)      ~ 1 | country_iso_alpha3:window,
      log(nbdisp) ~ 1 | country_iso_alpha3:window
    ),
    data_par    = endpoints,
    endpoints   = endpoints,
    curve       = "logistic",
    distr       = "nbinom",
    na_action   = c("exclude", "fail"),
    init        = init
  )
  save(object_m, file = "object_m.RData")

  pdf("world_incidence_m.pdf", width = 8, height = 4, onefile = TRUE)
  plot(object,
    type = "interval",
    show_tdoubling = TRUE,
    xlim = c(as.Date("2020-01-01"), Sys.Date()),
    log = TRUE,
    sub = country
  )
  dev.off()
}

sink(type = "output")
sink(type = "message")
