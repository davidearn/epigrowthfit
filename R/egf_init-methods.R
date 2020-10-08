#' Methods for class "egf_init"
#'
#' @description
#' Methods for "egf_init" objects returned by [egf_init()].
#'
#' @param x,object An "egf_init" object.
#' @param log A logical scalar. If `TRUE`, then parameter values are
#'   log-transformed.
#' @param time A numeric vector listing (increasing) time points in days
#'   since `object$date[1]`. Missing values are not tolerated. Must have
#'   length 2 or greater.
#' @param nsim A positive integer specifying a number of simulations.
#' @param seed An integer specifying a seed for RNG, otherwise `NULL`.
#' @param inc One of `"interval"` and `"cumulative"`,
#'   indicating a type of incidence to plot.
#' @param tol A non-negative number used only if `inc = "interval"`.
#'   `x$cases[i]` is plotted in light red if `diff(x$time)[i] < (1-tol)*m`,
#'   in dark red if `diff(x$time)[i] > (1+tol)*m`, and in grey otherwise,
#'   where `m = median(diff(x$time))`. Assign `Inf` to ensure that everything
#'   is grey.
#' @param ... Unused optional arguments.
#'
#' @return
#' The `print` method returns `x` invisibly.
#'
#' The `coef` method returns `object$theta0` if `log = FALSE`
#' and `object$log_theta0` if `log = TRUE`.
#'
#' The `predict` method returns a list with numeric elements:
#'
#' \describe{
#'   \item{`time`}{Matches argument.}
#'   \item{`cum_inc`}{Expected cumulative incidence at time points `time`,
#'     conditional on parameter values `object$theta0`.
#'     Equal to `object$cum_inc(time)`.
#'   }
#'   \item{`int_inc`}{Expected interval incidence given `time` as interval
#'     endpoints, conditional on parameter values `object$theta0`.
#'     Equal to `object$int_inc(time)`, which evaluates to `diff(cum_inc)`.
#'   }
#' }
#'
#' The `simulate` method returns a list with numeric elements:
#'
#' \describe{
#'   \item{`time`}{Matches argument.}
#'   \item{`cum_inc`}{A matrix with `length(time)` rows and `nsim` columns,
#'     with each column a cumulative incidence curve implied by a simulated
#'     interval incidence curve (see `int_inc` below). Column `i` is computed
#'     as `x0 + cumsum(c(0, int_inc[, i]))`
#'     where `x0 = object$cum_inc(time[1])` is the predicted value of
#'     cumulative incidence at `time[1]`, conditional on parameter values
#'     `object$theta0`.
#'   }
#'   \item{`int_inc`}{A matrix with `length(time)-1` rows and `nsim` columns,
#'     with each column a simulated interval incidence curve. Simulations are
#'     performed by adding observation error to `predict(object, time)$int_inc`,
#'     using the observation model specified by `object$distr`.
#'   }
#' }
#'
#' The `plot` method returns `NULL` invisibly.
#'
#' @details
#' ## Plot elements
#'
#' The bottom axis measures the number of days since `x$date[1]`.
#' The left axis measures interval or cumulative incidence (depending
#' on `inc`) on a log scale. Zeros are plotted as if they were 10^-0.2,
#' and are therefore distinguished from nonzero counts, which are always
#' at least 1.
#'
#' Observed data, specified by `x$time` and either `x$cases` or
#' `cumsum(x$cases)`, are plotted as filled points. `cases[i]` gives
#' the number of cases observed between `time[i]` and `time[i+1]`,
#' and `cumsum(cases)[i]` the number observed between `time[1]` and
#' `time[i+1]`. Both are plotted at `time[i+1]`.
#'
#' The left and right endpoints of the fitting window, specified by indices
#' `x$first` and `x$last`, are displayed as vertical lines at `time[first+1]`
#' and `time[last+1]` (adding one since `first` and `last` index `cases`,
#' and `length(time) = length(cases)+1`).
#'
#' The incidence curve predicted by initial parameter estimates `x$theta0`
#' is plotted as a teal line on grid points
#' `wgrid = seq(time[first+1], time[last+1], by = m)`,
#' where `m = median(diff(time))` for interval incidence (ensuring that
#' the curve has the correct scale; see below) and `m = 1` for cumulative
#' incidence. The predicted curve is obtained with `predict(x, wgrid)`.
#' The initial parameter estimates are printed at the bottom of the right
#' margin.
#'
#' Careful interpretation of observed interval incidence is required
#' if the plotted time series is not equally spaced, as `cases` roughly
#' scales with `diff(time)`. That is, certain observations may vary from
#' the predicted curve not due to chance, but because they represent a
#' count over fewer or more days than the typical observation interval,
#' namely `median(diff(time))`. Observations for which `diff(time)`
#' differs sufficently from `median(diff(time))` are highlighted according
#' to argument `tol`.
#'
#' @seealso [egf_init()]
#'
#' @name egf_init-methods
NULL

#' @rdname egf_init-methods
#' @export
print.egf_init <- function(x, ...) {
  if (!inherits(x, "egf_init")) {
    stop("`x` must be an \"egf_init\" object.")
  }

  cstr <- switch(x$curve,
    exponential = "an exponential model",
    logistic    = "a logistic model",
    richards    = "a Richards model"
  )
  bstr <- if (x$include_baseline) "with a linear baseline and" else "with"
  dstr <- switch(x$distr,
    poisson = "Poisson-distributed observations.",
    nbinom  = "negative binomial observations."
  )
  f <- x$first
  l <- x$last
  cat("Pass this \"egf_init\" object to `egf()` to fit", cstr, "\n")
  cat(bstr, dstr, "\n")
  cat("\n")
  cat("Fitting window:\n")
  cat("\n")
  cat("index   ", f, ":", l, "\n", sep = "")
  cat(" date   [", as.character(x$date[f+1]), ", ", as.character(x$date[l+1]), "]\n", sep = "")
  cat("cases   ", sum(x$cases[f:l]), " of ", sum(x$cases), "\n", sep = "")
  cat("\n")
  cat("Initial parameter estimates:\n")
  cat("\n")
  print(x$theta0)
  invisible(x)
}

#' @rdname egf_init-methods
#' @export
coef.egf_init <- function(object, log = FALSE, ...) {
  if (!inherits(object, "egf_init")) {
    stop("`object` must be an \"egf_init\" object.")
  }
  if (!is.logical(log) || length(log) != 1 || is.na(log)) {
    stop("`log` must be `TRUE` or `FALSE`.")
  }

  if (log) {
    object$log_theta0
  } else {
    object$theta0
  }
}

#' @rdname egf_init-methods
#' @export
predict.egf_init <- function(object, time, ...) {
  if (!inherits(object, "egf_init")) {
    stop("`object` must be an \"egf_init\" object.")
  }
  if (missing(time)) {
    stop("Missing argument `time`.")
  } else if (!is.numeric(time) || length(time) < 2) {
    stop("`time` must be numeric and have length 2 or greater.")
  } else if (anyNA(time)) {
    stop("`time` must not have missing values.")
  } else if (!all(diff(time) > 0)) {
    stop("`time` must be increasing.")
  }

  list(
    time = time,
    cum_inc = object$cum_inc(time),
    int_inc = object$int_inc(time)
  )
}

#' @rdname egf_init-methods
#' @export
#' @importFrom stats rpois rnbinom
simulate.egf_init <- function(object, nsim = 1, seed = NULL, time, ...) {
  if (!inherits(object, "egf_init")) {
    stop("`object` must be an \"egf_init\" object.")
  }
  if (missing(time)) {
    stop("Missing argument `time`.")
  } else if (!is.numeric(time) || length(time) < 2) {
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

  x <- object$cum_inc(time)
  dx <- diff(x)
  if (object$distr == "pois") {
    set.seed(seed)
    sim <- replicate(nsim, rpois(length(dx), dx))
  } else if (object$distr == "nbinom") {
    k <- object$theta0[["nbdisp"]]
    set.seed(seed)
    sim <- replicate(nsim, rnbinom(length(dx), mu = dx, size = k))
  }

  list(
    time = time,
    cum_inc = x[1] + rbind(0, apply(sim, 2, cumsum)),
    int_inc = sim
  )
}

#' @rdname egf_init-methods
#' @export
#' @import graphics
#' @importFrom stats median
plot.egf_init <- function(x, inc = "interval", tol = 0.025, ...) {
  if (!inherits(x, "egf_init")) {
    stop("`x` must be an \"egf_init\" object.")
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

  op <- par(mar = c(5, 4, 4, 8) + 0.1, las = 1, mgp = c(3, 0.7, 0))
  data <- data.frame(
    time = x$time[-1],
    int_inc = x$cases,
    cum_inc = cumsum(x$cases)
  )
  data$dt <- diff(x$time)
  m <- median(data$dt)
  dt_min <- (1 - tol) * m
  dt_max <- (1 + tol) * m
  bg_enum <- 1 + 1 * (data$dt < dt_min) + 2 * (data$dt > dt_max)
  data$bg <- c("#DDDDDD", "#CC6677", "#882255")[bg_enum]
  data$col <- c("#BBBBBB", "#CC6677", "#882255")[bg_enum]
  xlab <- paste("days since", as.character(x$date[1]))
  ylab <- "cases"
  f <- x$first
  l <- x$last
  ## Interval incidence
  if (inc == "interval") {
    data_z <- data[data$int_inc == 0, ]
    data_nz <- data[data$int_inc > 0, ]
    wgrid <- seq(x$time[f+1], x$time[l+1], by = m)
    pred <- predict(x, wgrid)
    ylim <- c(10^-0.2, max(c(data$int_inc, pred$int_inc)) * 10^0.2)
    yax_at <- 10^(0:max(floor(log10(c(data$int_inc, pred$int_inc)))))
    yax_labels <- parse(text = paste0("10^", log10(yax_at)))
    plot(int_inc ~ time, data = data_nz, ylim = ylim, yaxs = "i", yaxt = "n",
         log = "y", pch = 21, bg = data_nz$bg, col = data_nz$col,
         xlab = xlab, ylab = ylab)
    points(data_z$time, rep(ylim[1], nrow(data_z)), xpd = NA,
           pch = 21, bg = data_z$bg, col = data_z$col)
    lines(pred$time[-1], pred$int_inc, lwd = 3, col = "#44AA99")
    axis(side = 2, at = yax_at, labels = yax_labels)
  ## Cumulative incidence
  } else if (inc == "cumulative") {
    data_z <- data[data$cum_inc == 0, ]
    data_nz <- data[data$cum_inc > 0, ]
    wgrid <- seq(x$time[f+1], x$time[l+1], by = 1)
    pred <- predict(x, wgrid)
    ylim <- c(10^-0.2, max(c(data$cum_inc, pred$cum_inc)) * 10^0.2)
    yax_at <- 10^(0:max(floor(log10(c(data$cum_inc, pred$cum_inc)))))
    yax_labels <- parse(text = paste0("10^", log10(yax_at)))
    plot(cum_inc ~ time, data = data_nz, ylim = ylim, yaxs = "i", yaxt = "n",
         log = "y", pch = 21, bg = "#DDDDDD", col = "#BBBBBB",
         xlab = xlab, ylab = ylab)
    points(data_z$time, rep(ylim[1], nrow(data_z)), xpd = NA,
           pch = 21, bg = "#DDDDDD", col = "#BBBBBB")
    lines(wgrid, pred$cum_inc, lwd = 3, col = "#44AA99")
    axis(side = 2, at = yax_at, labels = yax_labels)
  }
  ## Fitting window
  abline(v = x$time[c(f,l)+1], lty = 2, col = "#555555")
  ## Initial parameter estimates
  pstr1 <- paste0(names(x$theta0), " = ")
  pstr2 <- round(x$theta0, digits = 4)
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
  if (inc == "interval") {
    legend(x = lx, y = ly, xpd = NA, bty = "n", cex = 0.7, seg.len = 1,
           pch = c(21, 21, 21, NA),
           pt.bg = c("#DDDDDD", "#CC6677", "#882255", NA),
           lty = c(NA, NA, NA, 1),
           lwd = c(NA, NA, NA, 3),
           col = c("#BBBBBB", "#CC6677", "#882255", "#44AA99"),
           legend = c("obs, dt ~ median",
                      "obs, dt < median",
                      "obs, dt > median",
                      "pred"))
  } else if (inc == "cumulative") {
    legend(x = lx, y = ly, xpd = NA, bty = "n", cex = 0.7, seg.len = 1,
           pch = c(21, NA),
           pt.bg = c("#DDDDDD", NA),
           lty = c(NA, 1),
           lwd = c(NA, 3),
           col = c("#BBBBBB", "#44AA99"),
           legend = c("obs", "pred"))
  }
  ## Title
  cstr <- x$curve
  substr(cstr, 1, 1) <- toupper(substr(cstr, 1, 1))
  title(main = paste(cstr, "model of", inc, "incidence\n(initial guess)"),
        cex.main = 0.9)
  par(op)
  invisible(NULL)
}
