#' Methods for class "egf"
#'
#' @description
#' Methods for "egf" objects returned by [egf()].
#'
#' @param x,object An "egf" object.
#' @param log A logical scalar. If `TRUE`, then parameter values are
#'   log-transformed.
#' @param time A numeric vector listing increasing time points in days
#'   since `object$init$date[1]`. Missing values are not tolerated.
#'   The simulate method requires `length(time) >= 2`.
#' @param nsim A positive integer specifying a number of simulations.
#' @param seed An integer specifying a seed for RNG, otherwise `NULL`.
#' @param inc One of `"interval"` and `"cumulative"`,
#'   indicating a type of incidence to plot.
#' @param tol A non-negative number used only if `inc = "interval"`.
#'   `x$init$cases[i]` is highlighted if
#'   `diff(x$init$time)[i] < (1-tol)*m` or
#'   `diff(x$init$time)[i] > (1+tol)*m`, where
#'   `m = median(diff(x$init$time))`.
#'   Assign 0 to ensure that all deviations from `m` are highlighted.
#'   Assign `Inf` to disable highlighting.
#' @param ... Unused optional arguments.
#'
#' @return
#' The `print` method returns `x` invisibly.
#'
#' The `coef` method returns `object$theta_hat` if `log = FALSE`
#' and `object$log_theta_hat` if `log = TRUE`.
#'
#' The `predict` method returns a list with numeric elements:
#'
#' \describe{
#'   \item{`time`}{Matches argument.}
#'   \item{`cum_inc`}{Expected cumulative incidence at time points
#'     `time`, conditional on parameter vector `object$theta_hat`.
#'     Equal to `object$cum_inc(time)`.
#'   }
#'   \item{`int_inc`}{Expected interval incidence given interval
#'     endpoints `time`, conditional on parameter vector
#'     `object$theta_hat`. Equal to `diff(object$cum_inc(time))`
#'     if `length(time) >= 2` and omitted otherwise.
#'   }
#' }
#'
#' The `simulate` method returns a list with numeric elements:
#'
#' \describe{
#'   \item{`time`}{Matches argument.}
#'   \item{`int_inc`}{A matrix with `length(time)-1` rows and `nsim`
#'     columns, such that `int_inc[i, j]` is the number of cases
#'     observed between `time[i]` and `time[i+1]` in simulation `j`
#'     of `nsim`. Row vector `int_inc[i, ]` is sampled from a
#'     Poisson or negative binomial distribution (depending on
#'     `object$init$distr`) with mean `predict(object, time)$int_inc[i]`.
#'     The negative binomial dispersion parameter is taken from
#'     `object$theta_hat[["nbdisp"]]`.
#'   }
#'   \item{`cum_inc`}{A matrix with `length(time)` rows and `nsim`
#'     columns, such that `cum_inc[i, j]` is the number of cases
#'     observed up to `time[i]` in simulation `j`. Column vector
#'     `cum_inc[, j]` is computed as `c0 + cumsum(c(0, int_inc[, j]))`,
#'     where `c0 = predict(object, time)$cum_inc[1]` is the expected
#'     value of cumulative incidence at `time[1]` conditional on
#'     parameter vector `object$theta_hat`.
#'   }
#' }
#'
#' The `plot` method returns `NULL` invisibly.
#'
#' @details
#' ## Plot elements
#'
#' The bottom axis measures the number of days since `x$init$date[1]`.
#' The left axis measures interval or cumulative incidence
#' (depending on `inc`) on a log scale. Zeros are plotted as if
#' they were 10^(-0.2), and are therefore distinguished from
#' nonzero counts, which are always at least 1.
#'
#' Observed data, specified by `x$init$time` and either `x$init$cases`
#' or `cumsum(x$init$cases)` (depending on `inc`), are plotted as
#' points. `cases[i]` gives the number of cases observed between
#' `time[i]` and `time[i+1]`, and `cumsum(cases)[i]` the number
#' observed between `time[1]` and `time[i+1]`. Both are plotted at
#' `time[i+1]`.
#'
#' The left and right endpoints of the fitting window, specified by
#' indices `x$init$first` and `x$init$last`, are displayed as vertical
#' lines at `time[first+1]` and `time[last+1]` (adding one since `first`
#' and `last` index `cases`, and `length(time) = length(cases)+1`).
#'
#' The incidence curve predicted by fitted parameter estimates
#' `x$theta_hat` is plotted as a teal line on grid points
#' `wgrid = seq(time[first+1], time[last+1], by = m)`,
#' where `m = 1` for cumulative incidence and
#' where `m = median(diff(time))` for interval incidence
#' (to ensure that the curve has the correct scale; see below).
#  The predicted curve is obtained with `predict(x, wgrid)`.
#  The fitted parameter estimates are printed at the bottom
#  of the right margin.
#'
#' Careful interpretation of observed interval incidence is required
#' if the plotted time series is not equally spaced, as `cases`
#' roughly scales with `diff(time)`. That is, certain observations
#' may vary from the predicted curve not due to chance, but because
#' they represent a count over fewer or more days than the typical
#' observation interval, namely `median(diff(time))`. Observations
#' for which `diff(time)` differs from `median(diff(time))` are
#' highlighted according to argument `tol`.
#'
#' @seealso [egf()]
#'
#' @name egf-methods
NULL

#' @rdname egf-methods
#' @export
print.egf <- function(x, ...) {
  if (!inherits(x, "egf")) {
    stop("`x` must be an \"egf\" object.")
  }

  cstr <- switch(x$init$curve,
    exponential = "an exponential model",
    logistic    = "a logistic model",
    richards    = "a Richards model"
  )
  bstr <- if (x$init$include_baseline) "with a linear baseline and" else "with"
  dstr <- switch(x$init$distr,
    poisson = "Poisson-distributed observations.",
    nbinom  = "negative binomial observations."
  )
  f <- x$init$first
  l <- x$init$last
  cat("This \"egf\" object fits", cstr, "\n")
  cat(bstr, dstr, "\n")
  cat("\n")
  cat("Fitting window:\n")
  cat("\n")
  cat("index   ", f, ":", l, "\n", sep = "")
  cat(" date   [", as.character(x$init$date[f+1]), ", ", as.character(x$init$date[l+1]), "]\n", sep = "")
  cat("cases   ", sum(x$init$cases[f:l]), " of ", sum(x$init$cases), "\n", sep = "")
  cat("\n")
  cat("Fitted parameter estimates:\n")
  cat("\n")
  print(x$theta_hat)
  invisible(x)
}

#' @rdname egf-methods
#' @export
coef.egf <- function(object, log = FALSE, ...) {
  if (!inherits(object, "egf")) {
    stop("`object` must be an \"egf\" object.")
  }
  if (!is.logical(log) || length(log) != 1 || is.na(log)) {
    stop("`log` must be `TRUE` or `FALSE`.")
  }

  if (log) {
    object$log_theta_hat
  } else {
    object$theta_hat
  }
}

#' @rdname egf-methods
#' @export
predict.egf <- function(object, time = object$init$time, ...) {
  if (!inherits(object, "egf")) {
    stop("`object` must be an \"egf\" object.")
  }
  if (!is.numeric(time) || length(time) == 0) {
    stop("`time` must be numeric and have nonzero length.")
  } else if (anyNA(time)) {
    stop("`time` must not have missing values.")
  } else if (!all(diff(time) > 0)) {
    stop("`time` must be increasing.")
  }

  out <- list(time = time, cum_inc = object$cum_inc(time))
  if (length(time) > 1) {
    out$int_inc = diff(out$cum_inc)
  }
  out
}

#' @rdname egf-methods
#' @export
#' @importFrom stats rpois rnbinom
simulate.egf <- function(object, nsim = 1, seed = NULL,
                         time = object$init$time, ...) {
  if (!inherits(object, "egf")) {
    stop("`object` must be an \"egf\" object.")
  }
  if (!is.numeric(time) || length(time) < 2) {
    stop("`time` must be numeric and have length 2 or greater.")
  } else if (anyNA(time)) {
    stop("`time` must not have missing values.")
  } else if (!all(diff(time) > 0)) {
    stop("`time` must be increasing.")
  }
  if (!is.numeric(nsim) || length(nsim) != 1 || !isTRUE(nsim >= 1)) {
    stop("`nsim` must be a positive integer.")
  }
  if (!is.null(seed)) {
    if (!is.null(numeric) && (!is.numeric(seed) || !is.finite(seed[1]))) {
      stop("`seed` must be `NULL` or an integer.")
    }
  }

  cum_inc <- object$cum_inc(time)
  int_inc <- diff(cum_inc)
  if (object$init$distr == "pois") {
    set.seed(seed)
    sim <- replicate(nsim, rpois(int_inc, lambda = int_inc))
  } else if (object$init$distr == "nbinom") {
    k <- object$theta_hat[["nbdisp"]]
    set.seed(seed)
    sim <- replicate(nsim, rnbinom(int_inc, mu = int_inc, size = k))
  }

  list(
    time = time,
    int_inc = sim,
    cum_inc = cum_inc[1] + apply(rbind(0, sim), 2, cumsum)
  )
}

#' @rdname egf-methods
#' @export
#' @import graphics
#' @importFrom stats median
plot.egf <- function(x, inc = "cumulative", tol = 0.025, ...) {
  if (!inherits(x, "egf")) {
    stop("`x` must be an \"egf\" object.")
  }
  if (!is.character(inc) || length(inc) != 1 ||
        !inc %in% c("interval", "cumulative")) {
    stop("`inc` must be \"interval\" or \"cumulative\".")
  }
  if (inc == "interval") {
    if (!is.numeric(tol) || length(tol) != 1 || !isTRUE(tol >= 0)) {
      stop("`tol` must be a non-negative number.")
    }
  }

  data <- data.frame(
    time = x$init$time[-1],
    int_inc = x$init$cases,
    cum_inc = cumsum(x$init$cases)
  )
  data$dt <- diff(x$init$time)
  m <- median(data$dt)
  dt_min <- (1 - tol) * m
  dt_max <- (1 + tol) * m
  dt_enum <- 1 + 1 * (data$dt < dt_min) + 2 * (data$dt > dt_max)
  data$bg <- c("#DDDDDD", "#FFFFFF", "#882255")[dt_enum]
  data$col <- c("#BBBBBB", "#882255", "#882255")[dt_enum]
  f <- x$init$first
  l <- x$init$last
  wgrid <- seq(x$init$time[f+1], x$init$time[l+1], by = 1)
  pred <- predict(x, wgrid)
  xlim <- c(0, max(x$init$time) * 1.04)
  xlab <- paste("days since", as.character(x$init$date[1]))
  ylab <- "cases"

  op <- par(mar = c(5, 4, 4, 8) + 0.1, las = 1, mgp = c(3, 0.7, 0))

  ## Cumulative incidence
  if (inc == "cumulative") {
    ylim <- c(10^-0.2, max(c(data$cum_inc, pred$cum_inc)) * 10^0.2)
    yax_at <- 10^(0:max(floor(log10(c(data$cum_inc, pred$cum_inc)))))
    yax_labels <- parse(text = paste0("10^", log10(yax_at)))
    data$cum_inc[data$cum_inc == 0] <- ylim[1]

    plot.new()
    plot.window(xlim = xlim, ylim = ylim,
                xaxs = "i", yaxs = "i", log = "y")
    axis(side = 1)
    axis(side = 2, at = yax_at, labels = yax_labels)
    box(bty = "l")
    points(cum_inc ~ time, data = data, xpd = NA,
           pch = 21, bg = "#DDDDDD", col = "#BBBBBB")
    lines(pred$time, pred$cum_inc, lwd = 3, col = "#44AA99")

    ## Interval incidence
  } else if (inc == "interval") {
    ylim <- c(10^-0.2, max(c(data$int_inc, pred$int_inc)) * 10^0.2)
    yax_at <- 10^(0:max(floor(log10(c(data$int_inc, pred$int_inc)))))
    yax_labels <- parse(text = paste0("10^", log10(yax_at)))
    data$int_inc[data$int_inc == 0] <- ylim[1]

    plot.new()
    plot.window(xlim = xlim, ylim = ylim,
                xaxs = "i", yaxs = "i", log = "y")
    axis(side = 1)
    axis(side = 2, at = yax_at, labels = yax_labels)
    box(bty = "l")
    points(int_inc ~ time, data = data, xpd = NA,
           pch = 21, bg = data$bg, col = data$col)
    lines(pred$time[-1], pred$int_inc, lwd = 3, col = "#44AA99")
  }

  ## Fitting window
  abline(v = x$init$time[c(f,l)+1], lty = 2, col = "#555555")

  ## Titles
  title(xlab = xlab)
  title(ylab = ylab)
  cstr <- x$init$curve
  substr(cstr, 1, 1) <- toupper(substr(cstr, 1, 1))
  title(main = paste(cstr, "model of", inc, "incidence\n(fitted)"),
        cex.main = 0.9)

  ## Initial parameter estimates
  pstr1 <- paste0(names(x$theta_hat), " = ")
  pstr2 <- round(x$theta_hat, digits = 4)
  if ("K" %in% names(x$theta_hat)) {
    pstr2[["K"]] <- round(pstr2[["K"]])
  }
  px <- par("usr")[2] + 0.02 * diff(par("usr")[1:2])
  px <- px + max(strwidth(pstr1, cex = 0.7))
  py <- 10^(par("usr")[3] + 0.02 * diff(par("usr")[3:4]))
  text(px, py, paste(pstr1, collapse = "\n"),
       adj = c(1, 0), xpd = NA, cex = 0.7)
  text(px, py, paste(pstr2, collapse = "\n"),
       adj = c(0, 0), xpd = NA, cex = 0.7)

  ## Legend
  lx <- par("usr")[2] + 0.02 * diff(par("usr")[1:2])
  ly <- 10^(par("usr")[4] - 0.02 * diff(par("usr")[3:4]))
  if (inc == "cumulative") {
    legend(x = lx, y = ly,
           xpd = NA, bty = "n", cex = 0.7, seg.len = 1,
           pch = c(21, NA),
           pt.bg = c("#DDDDDD", NA),
           lty = c(NA, 1),
           lwd = c(NA, 3),
           col = c("#BBBBBB", "#44AA99"),
           legend = c("obs", "pred")
    )
  } else if (inc == "interval") {
    legend(x = lx, y = ly,
           xpd = NA, bty = "n", cex = 0.7, seg.len = 1,
           pch = c(21, 21, 21, NA),
           pt.bg = c("#DDDDDD", "#FFFFFF", "#882255", NA),
           lty = c(NA, NA, NA, 1),
           lwd = c(NA, NA, NA, 3),
           col = c("#BBBBBB", "#882255", "#882255", "#44AA99"),
           legend = c("obs, dt ~ median",
                      "obs, dt < median",
                      "obs, dt > median",
                      "pred")
    )
  }

  par(op)
  invisible(NULL)
}
