f <- function(formula, data, spar) {
  v <- all.vars(formula)
  date <- data[[v[2L]]]
  cases <- data[[v[1L]]]
  data <- na.omit(data.frame(
    date = date,
    time = as.integer(date - date[1L]),
    difftime = c(NA_integer_, as.integer(diff(date))),
    cases = cases,
    log1p_cases = log10(1L + cases)
  ))

  spar <- 0.6
  ss <- smooth.spline(
    x = data$time,
    y = data$log1p_cases,
    w = data$difftime,
    spar = spar
  )
  dss <- predict(ss, x = data$time, deriv = 1L)

  op <- par(mfrow = c(2L, 1L), mar = c(2, 4, 0.5, 0.5) + 0.1,
            cex.axis = 0.8, cex.lab = 0.8, mgp = c(3, 0.7, 0), las = 1)
  plot(log1p_cases ~ time, data = data)
  lines(ss, lwd = 2, col = "red")
  plot(log(y) ~ x, data = dss, type = "l", lwd = 2, col = "blue",
       ylab = "log(d_log1p_cases)")
  par(op)
  invisible(NULL)
}

f(new_confirmed ~ date, data = ontario, spar = 0.6)
f(new_confirmed ~ date, data = quebec,  spar = 0.6)



