library("epigrowthfit")
load("../world/data/world.RData")
load("../world/data/endpoints.RData")
options(contrasts = c("contr.sum", "contr.poly"), warn = 1)
f <- function(x) { # dummy to contr.sum
  m <- mean(x)
  c(m, x[-length(x)] - m)
}

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
  Y_init <- rbind(Y_init, matrix(object_1ts$report$Y_as_vector, ncol = 4L))
}

## Initial parameter vector for mixed effects model
m <- colMeans(Y_init)
beta <- c(f(Y_init[, 1L]), m[2:4])
o <- order(endpoints$window, endpoints$country_iso_alpha3)
b <- t(Y_init[o, 2:4]) - m[2:4]
log_sd_b <- log(apply(matrix(b, nrow = 3L), 1L, sd))
b <- b / exp(log_sd_b)
theta <- c(log_sd_b, rep_len(0, 3L))
init <- c(beta, b, theta)

## Mixed effects model fit
outfile <- file("egf.Rout", open = "wt")
sink(outfile, type = "output")
sink(outfile, type = "message")
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
  trace       = 3L,
  init        = init
)
sink(type = "message")
sink(type = "output")

object$tmb_out$fn(object$init)
debugonce(optimHess)
sdreport(object$tmb_out)

options(warn = 2L, error = recover)
options(warn = 1L, error = NULL)
with(object$tmb_out,
  nlminb(env$last.par.best[-env$random], fn, gr, control = list(trace = 1L))
)
zz <- with(object$tmb_out,
  optim(env$last.par.best[-env$random], fn, gr, method = "BFGS", control = list(trace = 1L))
)
sdreport(object$tmb_out, par.fixed = zz$par)

oh <- with(object$tmb_out,
  optimHess(env$last.par.best[-env$random], fn, gr, control = list(trace = 100L))
)
summary(log10(abs(object$best)))
oh_parscale <- with(object$tmb_out,
  optimHess(env$last.par.best[-env$random], fn, gr,
            control = list(trace = 100L,
                           parscale = abs(env$last.par.best[-env$random])))
)
fn_par <- gr_par <- NULL
fn2 <- function(x, ...) {
  fn_par <<- x
  print(x)
  res <- object$tmb_out$fn(x, ...)
  cat("fn = ", res, "\n")
  res
}
gr2 <- function(x, ...) {
  gr_par <<- x
  print(x)
  res <- object$tmb_out$gr(x, ...)
  cat("gr = ", res, "\n")
  res
}


oh_parscale2 <- with(object$tmb_out,
  optimHess(env$last.par.best[-env$random], fn2, gr2,
            control = list(trace = 100L,
                           parscale = abs(env$last.par.best[-env$random])))
)
debugonce(object$tmb_out$gr)
object$tmb_out$env$inner.method <- "BFGS"
object$tmb_out$gr(gr_par)

# numDeriv::jacobian/hessian

save(object, file = "object.RData")
