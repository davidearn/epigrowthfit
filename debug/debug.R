library("epigrowthfit")
load("../world/data/world.RData")
load("../world/data/endpoints.RData")
options(contrasts = c("contr.sum", "contr.poly"), warn = 1)
f <- function(x) { # dummy to contr.sum
  m <- mean(x)
  c(m, x[-length(x)] - m)
}

## Want to fit mixed effects model
##
## log(r)      ~ 1 + I(country:window)
## log(tinfl)  ~ 1 + (1 | country:window)
## log(K)      ~ 1 + (1 | country:window)
## log(nbdisp) ~ 1 + (1 | country:window)
##
## Need to start with a decent initial parameter vector
## to mitigate convergence issues. Idea will be to use
## (appropriate linear combinations of) fitted values from
## the fixed effects model
##
## log(r)      ~ 1 + I(country:window)
## log(tinfl)  ~ 1 + I(country:window)
## log(K)      ~ 1 c+ I(country:window)
## log(nbdisp) ~ 1 + I(country:window)
##
## However, the fixed effects model has convergence issues
## of its own when the number of fitting windows (levels of
## country:window) is large. A fast solution is to fit each
## time series separately and combine the resulting fitted
## values.

## Matrix of fitted values of log(r), log(tinfl), log(K), log(nbdisp)
## from fixed effects model
Y_init <- matrix(numeric(0L), ncol = 4L)
endpoints_split <- split(endpoints, endpoints$country_iso_alpha3)
cia3 <- names(endpoints_split)
for (i in seq_along(endpoints_split)) {
  print(cia3[i])
  object_1ts <- egf(
    formula_ts  = cases_new ~ Date | country_iso_alpha3,
    data_ts     = world,
    formula_par = if (nrow(endpoints_split[[i]]) > 1L) ~window else ~1,
    data_par    = endpoints_split[[i]],
    endpoints   = endpoints_split[[i]],
    curve       = "logistic",
    distr       = "nbinom",
    na_action   = c("exclude", "fail"),
    debug       = TRUE
  )
  init <- unlist(lapply(object_1ts$Y_init, f), use.names = FALSE)
  object_1ts <- update(object_1ts, debug = FALSE, init = init)
  Y_init <- rbind(Y_init, matrix(fitted(object_1ts)$estimate, ncol = 4L))
}

## Initial parameter vector for mixed effects model
m <- colMeans(Y_init)
beta <- c(f(Y_init[, 1L]), m[2:4])
b <- Y_init[, 2:4] - rep(m[2:4], each = nrow(Y_init))
log_sd_b <- log(apply(matrix(b, ncol = 3L), 2L, sd))
#b <- rep_len(0, 3L * nrow(Y_init))
#log_sd_b <- rep_len(1, 3L)
init <- c(beta, b, log_sd_b)

## Mixed effects model fit
debugonce(epigrowthfit:::optim_tmb_out)
object <- egf(
  formula_ts  = cases_new ~ Date | country_iso_alpha3,
  data_ts     = world,
  formula_par = list(
    log(r)      ~ I(country_iso_alpha3:window),
    log(tinfl)  ~ 1 | country_iso_alpha3:window,
    log(K)      ~ 1 | country_iso_alpha3:window,
    log(nbdisp) ~ 1 | country_iso_alpha3:window
  ),
  data_par    = endpoints,
  endpoints   = endpoints,
  curve       = "logistic",
  distr       = "nbinom",
  na_action   = c("exclude", "fail"),
  trace       = 1L,
  init        = init
)
outfile <- file("debug.Rout", open = "wt")
sink(outfile, type = "output")
sink(outfile, type = "message")
tmb_out$fn(tmb_out$par)
sink(type = "output")
sink(type = "message")
Q




o <- do.call(order, unname(endpoints[c("window", "country_iso_alpha3")]))
b_nolex <- Y_init[o, 2:4] - rep(m[2:4], each = nrow(Y_init))
init_nolex <- c(beta, b_nolex, log_sd_b)

outfile_nolex <- file("debug_nolex.Rout", open = "wt")
sink(outfile_nolex, type = "output")
sink(outfile_nolex, type = "message")
object <- egf(
  formula_ts  = cases_new ~ Date | country_iso_alpha3,
  data_ts     = world,
  formula_par = list(
    log(r)      ~ I(country_iso_alpha3:window),
    log(tinfl)  ~ 1 | country_iso_alpha3:window,
    log(K)      ~ 1 | country_iso_alpha3:window,
    log(nbdisp) ~ 1 | country_iso_alpha3:window
  ),
  data_par    = endpoints,
  endpoints   = endpoints,
  curve       = "logistic",
  distr       = "nbinom",
  na_action   = c("exclude", "fail"),
  trace       = 1L,
  init        = init_nolex
)
sink(type = "output")
sink(type = "message")
