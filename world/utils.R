do_plot_ <- function(formula, data, spar = NULL, endpoints = NULL, ...) {
  op <- par(mar = c(2, 3, 1, 1), mgp = c(2, 0.7, 0), las = 1)
  on.exit(par(op))
  xat <- seq(as.Date("2020-01-01"), Sys.Date(), by = "1 month")
  plot(formula, data = data, xlab = "", xaxt = "n", ...)
  abline(v = xat, lty = 3, col = "grey80")
  axis(side = 1, at = xat, labels = months(xat, TRUE), gap.axis = 0)
  if (!is.null(spar)) {
    lhs <- eval(formula[[2L]], envir = data, enclos = parent.frame())
    rhs <- eval(formula[[3L]], envir = data, enclos = parent.frame())
    dd <- data.frame(x = rhs - rhs[1], y = lhs, w = 1 / c(NA, diff(rhs)))
    dd <- dd[complete.cases(dd), , drop = FALSE]
    ss <- smooth.spline(x = dd$x, y = dd$y, w = dd$w, spar = spar)
    lines(y ~ I(rhs[1] + x), data = ss[c("x", "y")], lwd = 2, col = "red")
  }
  if (!is.null(endpoints)) {
    abline(v = as.Date(unlist(endpoints)), lty = 3)
  }
  invisible(NULL)
}
do_plot <- function(country_iso_alpha3, spar = NULL, endpoints = NULL, ...) {
  formula <- log1p(cases / dt) ~ date
  data <- world7_split[[country_iso_alpha3]]
  data$dt <- c(7, diff(data$date))
  dt7 <- data$dt == 7
  do_plot_(formula,
    data = data,
    spar = spar[2L],
    endpoints = endpoints,
    pch = 21,
    bg = c("blue", "grey80")[1L + dt7],
    main = sprintf("%s (weekly)", country_iso_alpha3),
    ...
  )
  if (!all(dt7)) {
    text(formula, data = data, subset = !dt7, labels = data$dt,
         pos = 3, offset = 0.3, col = "grey60", cex = 0.7)
  }
  data <- world_split[[country_iso_alpha3]]
  data$dt <- c(1, diff(data$date))
  dt1 <- data$dt == 1
  do_plot_(formula,
    data = data,
    spar = spar[1L],
    endpoints = endpoints,
    pch = 21,
    bg = c("blue", "grey80")[1L + dt1],
    main = sprintf("%s (daily)", country_iso_alpha3),
    ...
  )
  if (!all(dt1)) {
    text(formula, data = data, subset = !dt1, labels = data$dt,
         pos = 3, offset = 0.3, col = "grey60", cex = 0.7)
  }
  invisible(NULL)
}
