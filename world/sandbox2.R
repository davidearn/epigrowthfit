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

cia3 <- "SEN"
ep <- make_endpoints(
  list( # Senegal
    c("2020-03-05", "2020-03-25"),
    c("2020-04-11", "2020-04-30"),
    c("2020-11-21", "2020-12-11"),
    c("2020-12-13", "2021-01-14")
  )
)
{
d <- world[world$country_iso_alpha3 == cia3, , drop = FALSE]
object <- egf(
  formula_ts  = cases_new ~ Date,
  data_ts     = d,
  formula_par = if (nrow(ep) > 1L) ~window else ~1,
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
  #xlim = c("2020-02-01", "2020-08-01"),
  #xlim = c("2020-10-01", "2021-05-15"),
  #xlim = c("2021-01-01", "2021-06-01"),
  #ylim = c(50,2000)
)
print(object$optim_out$convergence)
}
fitted(object,se=TRUE)

