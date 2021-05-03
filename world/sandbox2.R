library("epigrowthfit")
load("data/world.RData")
load("data/endpoints.RData")
options(contrasts = c("contr.sum", "contr.poly"))

make_endpoints <- function(l, ...) {
  m <- matrix(unlist(l), ncol = 2L, byrow = TRUE,
              dimnames = list(NULL, c("start", "end")))
  d <- as.data.frame(m)
  d[] <- lapply(d, as.Date, ...)
  d <- d[order(d$start), , drop = FALSE]
  data.frame(d, window = gl(nrow(d), 1L))
}

cia3 <- "IND"
ep <- make_endpoints(
  list( # India
    c("2020-03-06", "2020-04-12"),
    c("2020-06-05", "2020-08-05"),
    c("2020-08-22", "2020-09-12"),
    c("2021-03-30", "2021-04-30")
  )
)
{
d <- world[world$country_iso_alpha3 == cia3, , drop = FALSE]
object <- egf(
  formula_ts  = cases_new ~ Date,
  data_ts     = d,
  #formula_par = ~1,
  formula_par = ~window,
  data_par    = ep,
  endpoints   = ep,
  curve       = "logistic",
  distr       = "nbinom",
  na_action   = c("exclude", "fail"),
  debug       = TRUE
)
## FIXME:
## Below construction of `init` assumes that
## match(levels(endpoints$country_iso_alpha3), levels(world$country_iso_alpha3), 0L)
## is increasing
f <- function(x) {
  m <- mean(x)
  c(m, x[-length(x)] - m)
}
init <- unlist(lapply(object$Y_init, f), use.names = FALSE)
object <- update(object, debug = FALSE, init = init)
plot(object,
  type = "interval",
  show_tdoubling = TRUE,
  log = TRUE,
  sub = cia3,
  xlim = c("2020-03-01", "2020-12-01"),
  #xlim = c("2020-08-01", "2021-05-01"),
  #xlim = c("2021-01-01", "2021-08-01"),
  ylim = c(100,500000)
)
print(object$optim_out$convergence)
}

