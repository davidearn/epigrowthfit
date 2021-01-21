f <- function(formula, data, spar, xlim = NULL) {
  v <- all.vars(formula)
  date <- data[[v[2L]]]
  cases <- data[[v[1L]]]
  data <- na.omit(data.frame(
    date = date,
    time = as.integer(date - date[1L]),
    difftime = c(NA_integer_, as.integer(diff(date))),
    cases = cases,
    log1p_cases = log1p(cases)
  ))

  ss <- smooth.spline(
    x = data$time,
    y = data$log1p_cases,
    w = data$difftime,
    spar = spar
  )
  dss <- as.data.frame(predict(ss, deriv = 1L))
  ddss <- as.data.frame(predict(ss, deriv = 2L))

  if (!is.null(xlim)) {
    data <- subset(data, time >= xlim[1L] & time <= xlim[2L])
    dss <- subset(dss, x >= xlim[1L] & x <= xlim[2L])
    ddss <- subset(ddss, x >= xlim[1L] & x <= xlim[2L])
  }


  op <- par(mfrow = c(3L, 1L), mar = c(2, 4, 0.5, 0.5) + 0.1,
            cex.axis = 0.8, cex.lab = 0.8, mgp = c(3, 0.7, 0), las = 1)
  plot(log1p_cases ~ time, data = data)
  lines(ss, lwd = 2, col = "red")
  plot(y ~ x, data = dss, type = "l", lwd = 2, col = "blue",
       ylab = "d1_log1p_cases")
  abline(h = 0)
  plot(y ~ x, data = ddss, type = "l", lwd = 2, col = "pink",
       ylab = "d2_log1p_cases")
  abline(h = 0)
  par(op)
  invisible(NULL)
}

f(new_confirmed ~ date, data = ontario, spar = 0.6, xlim = c(0,100))
f(new_confirmed ~ date, data = ontario, spar = 0.6, xlim = c(180,280))
f(new_confirmed ~ date, data = ontario, spar = 0.6, xlim = c(250,400))
f(new_confirmed ~ date, data = quebec,  spar = 0.7)
