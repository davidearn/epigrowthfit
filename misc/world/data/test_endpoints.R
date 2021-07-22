library("epigrowthfit")
load("../world.RData")
load("../endpoints.RData")
options(contrasts = c("contr.sum", "contr.poly"))
f <- function(x) { # dummy to zero-sum
  m <- mean(x)
  c(m, x[-length(x)] - m)
}
make_endpoints <- function(l, ...) {
  m <- matrix(unlist(l), ncol = 2L, byrow = TRUE,
              dimnames = list(NULL, c("start", "end")))
  d <- as.data.frame(m)
  d[] <- lapply(d, as.Date, ...)
  d <- d[order(d$start), , drop = FALSE]
  data.frame(d, window = gl(nrow(d), 1L))
}

cia3 <- "KOR"
ep <- make_endpoints(
  list( # South Korea
    c("2020-02-07", "2020-02-29"),
    c("2020-05-01", "2020-06-06"),
    c("2020-08-08", "2020-08-18"),
    c("2020-11-09", "2020-12-21"),
    c("2021-03-22", "2021-04-17")
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
init <- unlist(lapply(object$Y_init, f), use.names = FALSE)
object <- update(object, debug = FALSE, init = init)
plot(object,
  type = "interval",
  show_tdoubling = TRUE,
  log = TRUE,
  sub = cia3,
  #xlim = c("2020-02-01", "2020-10-01"),
  xlim = c("2020-08-01", "2021-06-01"),
  #xlim = c("2021-01-01", "2021-06-01"),
  #ylim = c(10, 1000)
)
# fitted(object, se = TRUE)
print(object$optim_out$convergence)
}


