library("epigrowthfit")
load("debug.RData")
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
## log(K)      ~ 1 + I(country:window)
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
for (d in split(endpoints, endpoints$country_iso_alpha3)) {
  object_1ts <- egf(
    formula_ts  = cases_new ~ Date | country_iso_alpha3,
    data_ts     = world,
    formula_par = if (nrow(d) > 1L) ~window else ~1,
    data_par    = d,
    endpoints   = d,
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
## Mixed effects model fit
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
  init        = c(beta, b, log_sd_b)
)
