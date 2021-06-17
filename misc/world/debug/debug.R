library("epigrowthfit")
load("../data/world.RData")
load("../data/endpoints.RData")
options(contrasts = c("contr.sum", "contr.poly"), warn = 1)
f <- function(x) { # dummy to contr.sum
  m <- mean(x)
  c(m, x[-length(x)] - m)
}

## Retrieve fitted values from within-country fixed effects models
## for initialization of across-countries mixed effects model
Y_init <- matrix(numeric(0L), ncol = 4L)
N <- c(table(endpoints$country_iso_alpha3))
cia3 <- names(N)
for (i in seq_along(N)) {
  print(cia3[i])
  object_1ts <- egf(
    formula = cases_new ~ Date | country_iso_alpha3,
    formula_par = if (N[[i]] > 1L) ~window else ~1,
    data = world,
    data_par = endpoints,
    subset_par = country_iso_alpha3 == cia3[i],
    na_action = "pass",
    endpoints = endpoints,
    do_fit = FALSE
  )
  init <- c(apply(object_1ts$Y_init, 2L, f))
  object_1ts <- update(object_1ts, do_fit = TRUE, init = init)
  Y_init <- rbind(Y_init, matrix(object_1ts$tmb_out$report(object_1ts$best)$Y_as_vector, ncol = 4L))
}


## Model without covariates

m <- colMeans(Y_init)
beta <- m
o <- order(endpoints$window, endpoints$country_iso_alpha3)
b <- t(Y_init[o, , drop = FALSE]) - m
log_sd_b <- log(apply(matrix(b, nrow = 4L), 1L, sd))
b <- b / exp(log_sd_b)
theta <- c(log_sd_b, rep_len(0, 6L))
init <- c(beta, b, theta)

outfile <- file("simple.Rout", open = "wt")
sink(outfile, type = "output")
sink(outfile, type = "message")
object <- egf(
  formula = cases_new ~ Date | country_iso_alpha3,
  formula_par = ~(1 | country_iso_alpha3:window),
  data = world,
  data_par = endpoints,
  na_action = "pass",
  endpoints = endpoints,
  init = init,
  se = TRUE,
  control = egf_control(
    trace = 2L,
    omp_num_threads = 4L,
    optimizer = egf_optimizer(
      f = optim,
      args = list(
        method = "BFGS"
      ),
      control = list(
        trace = 1L,
        maxit = 1000L
      )
    ),
    inner_optimizer = egf_inner_optimizer(
      f = newton,
      args = list(
        trace = 1L,
        maxit = 1000L
      )
    )
  )
)
sink(type = "message")
sink(type = "output")
save(object, file = "simple.RData")



library("tmbstan")
init_func <- function(m) {
  l <- split(m, sub("\\[[0-9]+\\]$", "", names(m)))
  function() {
    list(
      beta = rnorm(length(l$beta), l$beta),
      b = rnorm(l$b),
      theta = c(rnorm(4L, l$theta[1:4]), rnorm(l$theta[-(1:4)]))
    )
  }
}
zz <- tmbstan(object$tmb_out,
  init = object$best[object$nonrandom],
  #init = init_func(object$best),
  cores = 4L,
  laplace = TRUE,
  silent = FALSE
)


### Model with covariates

# tt <- substitute(list(
#   (1 | country_iso_alpha3:window),
#   log(mobility_retail_and_recreation),
#   log(mobility_grocery_and_pharmacy),
#   log(mobility_parks),
#   log(mobility_transit_stations),
#   log(mobility_workplaces),
#   log(mobility_residential),
#   qlogis(0.001 + npi_index_containment_health),
#   qlogis(0.001 + vaccinated_geq1d),
#   qlogis(econ_gini),
#   log(econ_gdp_per_capita),
#   log(weather_temperature),
#   log(weather_specific_humidity),
#   log(weather_shortwave_radiation),
#   log(weather_precipitation),
#   log(weather_wind_speed),
#   log(abs(latitude)),
#   days_since_2in1m
# ))
# plus <- function(x, y) call("+", x, y)
# sum_tt <- Reduce(plus, tt[-1L])
#
# cc <- complete.cases(endpoints[all.vars(tt)])
#
# m <- colMeans(Y_init[cc, , drop = FALSE])
# beta <- c(m[1L], rep_len(0, length(tt) - 2L), m[2:4])
# o <- order(endpoints$window[cc], endpoints$country_iso_alpha3[cc])
# b <- t(Y_init[which(cc)[o], , drop = FALSE]) - m
# log_sd_b <- log(apply(matrix(b, nrow = 4L), 1L, sd))
# b <- b / exp(log_sd_b)
# theta <- c(log_sd_b, rep_len(0, 6L))
# init <- c(beta, b, theta)
#
# outfile <- file("covariates.Rout", open = "wt")
# sink(outfile, type = "output")
# sink(outfile, type = "message")
# object_covariates <- update(object,
#   formula_par = list(
#     as.formula(call("~", quote(log(r)),      sum_tt)),
#     as.formula(call("~", quote(log(tinfl)),  tt[[2L]])),
#     as.formula(call("~", quote(log(K)),      tt[[2L]])),
#     as.formula(call("~", quote(log(nbdisp)), tt[[2L]]))
#   ),
#   subset_par = cc,
#   init = init
# )
# sink(type = "message")
# sink(type = "output")
