library("epigrowthfit")
load("../on_by_region.RData")
load("endpoints.RData")

window <- make_window(
  time = on_by_region$date,
  ts = on_by_region$region,
  endpoints = endpoints
)
wave <- make_wave(
  window = window,
  ts = on_by_region$region
)

object <- egf(
  formula_ts  = cases ~ date | region,
  formula_par = ~region + wave,
  data        = cbind(on_by_region, wave),
  window      = window,
  curve       = "logistic",
  obs_distr   = "nbinom"
)

ci <- confint(object, par = "log_r", link = FALSE, append = c(region, wave))
endpoints <- epigrowthfit:::get_window_endpoints(object)
endpoints[c("start", "end")] <- lapply(endpoints[c("start", "end")], `+`, attr(endpoints, "origin"))
on_r <- cbind(
  ci[c("region", "wave")],
  endpoints[c("start", "end")],
  ci[c("estimate", "lower", "upper")]
)
#save(on_td, file = "on_td.RData")


#pdf("on_td1.pdf", width = 6, height = 4, onefile = TRUE)
plot(ci, type = "bars", order = order(wave, region, estimate), label = wave:region)
#dev.off()

#pdf("on_td2.pdf", width = 6, height = 4, onefile = TRUE)
plot(ci, type = "boxes", label = region)
#dev.off()

# subset <- list(region = census_divisions)
# xlim <- c("2020-03-01", "2021-02-02")
#
# pdf("on_incidence.pdf", width = 8, height = 4, onefile = TRUE)
# plot(object,
#   type = "interval",
#   bands = TRUE,
#   subset = subset,
#   tol = Inf,
#   xlim = xlim,
#   ylab = "reported incidence",
#   main = "Daily reported COVID-19 incidence\n%region"
# )
# dev.off()
#
# pdf("on_rt1.pdf", width = 8, height = 4, onefile = TRUE)
# plot(object,
#   type = "rt1",
#   bands = TRUE,
#   subset = subset,
#   xlim = xlim,
#   ylim = c(0, log(2) / 2.5),
#   main = "Instantaneous exponential growth rate\n%region",
#   control = list(points = NULL, abline = NULL)
# )
# dev.off()
#
# pdf("on_rt2.pdf", width = 6, height = 4, onefile = TRUE)
# plot(object,
#   type = "rt2",
#   per_plot = 12,
#   bw_panels = 0,
#   subset = subset,
#   control = list(colorRamp = list(bias = 2))
# )
# dev.off()
