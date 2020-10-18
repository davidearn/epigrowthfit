#' Methods for class "egf_init"
#'
#' @description
#' Methods for "egf_init" objects returned by [egf_init()].
#'
#' @param x,object An "egf_init" object.
#' @param log A logical scalar. If `TRUE`, then parameter values are
#'   log-transformed.
#' @param time A numeric vector listing increasing time points in days
#'   since `object$date[1]`. Missing values are not tolerated.
#' @param inc One of `"cumulative"` and `"interval"`,
#'   indicating a type of incidence to plot.
#' @param tol A non-negative number used only if `inc = "interval"`.
#'   `x$cases[i]` is highlighted if
#'   `diff(x$time)[i] < (1-tol)*m`, or
#'   `diff(x$time)[i] > (1+tol)*m`, where
#'   `m = median(diff(x$time))`.
#'   Assign 0 to ensure that all deviations from `m` are highlighted.
#'   Assign `Inf` to disable highlighting.
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
#'   \item{`refdate`}{Matches `object$date[1]`.}
#'   \item{`cum_inc`}{Expected cumulative incidence at time points
#'     `time`, conditional on parameter vector `object$theta0`.
#'     Equal to `object$cum_inc(time)`.
#'   }
#'   \item{`int_inc`}{Expected interval incidence given interval
#'     endpoints `time`, conditional on parameter vector
#'     `object$theta0`. Equal to `diff(object$cum_inc(time))`
#'     if `length(time) >= 2` and omitted otherwise.
#'   }
#' }
#'
#' The `plot` method returns `NULL` invisibly.
#'
#' @details
#' ## Plot elements
#'
#' The bottom axis measures the number of days since `x$date[1]`.
#' The left axis measures interval or cumulative incidence
#' (depending on `inc`) on a log scale. Zeros are plotted as if
#' they were `10^(-0.2)`, and are therefore distinguished from
#' nonzero counts, which are always at least 1.
#'
#' Observed data, specified by `x$time` and either `x$cases` or
#' `cumsum(x$cases)` (depending on `inc`), are plotted as points.
#' `cases[i]` gives the number of cases observed between `time[i]`
#' and `time[i+1]`, and `cumsum(cases)[i]` the number observed
#' between `time[1]` and `time[i+1]`. Both are plotted at `time[i+1]`.
#'
#' The left and right endpoints of the fitting window, specified by
#' indices `x$first` and `x$last`, are displayed as vertical lines
#' at `time[first+1]` and `time[last+1]` (adding one since `first`
#' and `last` index `cases`, and `length(time) = length(cases)+1`).
#'
#' The incidence curve predicted by initial parameter estimates
#' `x$theta0` is plotted as a teal line on grid points
#' `wgrid = seq(time[first+1], time[last+1], by = m)`,
#' where `m = 1` for cumulative incidence and
#' where `m = median(diff(time))` for interval incidence
#' (to ensure that the curve has the correct scale; see below).
#  The predicted curve is obtained with `predict(x, wgrid)`.
#  The initial parameter estimates are printed at the bottom
#  of the right margin.
#'
#' Careful interpretation of observed interval incidence is required
#' if the plotted time series is not equally spaced, as `cases`
#' roughly scales with `diff(time)`. That is, certain observations
#' may vary from the predicted curve not due to chance, but because
#' they represent a count over fewer or more days than the typical
#' observation interval, namely `median(diff(time))`. Observations
#' for which `diff(time)` differs from `median(diff(time))` are
#' highlighted according to argument `tol` and labeled with the
#' value of `diff(time)`.
#'
#' @seealso [egf_init()]
#'
#' @name egf_init-methods
NULL

#' @rdname egf_init-methods
#' @export
print.egf_init <- function(x, ...) {
  cstr <- switch(x$curve,
    exponential = "an exponential model",
    logistic    = "a logistic model",
    richards    = "a Richards model"
  )
  bstr <- if (x$include_baseline) {
    "with a linear baseline and"
  } else {
    "with"
  }
  dstr <- switch(x$distr,
    poisson = "Poisson-distributed observations.",
    nbinom  = "negative binomial observations."
  )
  uvec <- c(r = "per day", thalf = "days", b = "per day")
  uvec <- uvec[names(uvec) %in% names(x$theta0)]
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
  cat("\n")
  cat("Units:\n")
  cat("\n")
  print(uvec, quote = FALSE)
  invisible(x)
}

#' @rdname egf_init-methods
#' @export
coef.egf_init <- function(object, log = FALSE, ...) {
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
predict.egf_init <- function(object, time = object$time, ...) {
  if (!is.numeric(time) || length(time) == 0) {
    stop("`time` must be numeric and have nonzero length.")
  } else if (anyNA(time)) {
    stop("`time` must not have missing values.")
  } else if (!all(diff(time) > 0)) {
    stop("`time` must be increasing.")
  }

  out <- list(
    time = time,
    refdate = object$date[1],
    cum_inc = object$cum_inc(time)
  )
  if (length(time) > 1) {
    out$int_inc = diff(out$cum_inc)
  }
  out
}

#' @rdname egf_init-methods
#' @export
#' @import graphics
#' @importFrom stats median
plot.egf_init <- function(x, inc = "cumulative", tol = 0, ...) {
  if (!is.character(inc) || length(inc) != 1 ||
      !inc %in% c("interval", "cumulative")) {
    stop("`inc` must be one of \"interval\", \"cumulative\".")
  }
  if (inc == "interval") {
    if (!is.numeric(tol) || length(tol) != 1 || !isTRUE(tol >= 0)) {
      stop("`tol` must be a non-negative number.")
    }
  }


  ### SET UP ###########################################################

  ## Optional graphical parameters
  dots <- list(...)

  ## Observed data
  data <- data.frame(
    time = x$time[-1],
    cum_inc = cumsum(x$cases),
    int_inc = x$cases
  )

  ## Predicted curve
  f <- x$first
  l <- x$last
  wgrid <- seq(x$time[f+1], x$time[l+1], by = 1)
  wpred <- predict(x, wgrid)[c("time", "cum_inc", "int_inc")]
  wpred$int_inc <- c(NA, wpred$int_inc)

  ## A way to avoid conditional `if (inc = ...) ... else ...`
  varname <- substr(inc, start = 1, stop = 3) # first three characters
  varname <- paste0(varname, "_inc")
  formula <- as.formula(paste(varname, "~ time"))

  ## Axis titles
  xlab <- paste("days since", as.character(x$date[1]))
  ylab <- paste(inc, "incidence")

  ## Axis limits (x)
  xmin <- 0
  xmax <- max(x$time, na.rm = TRUE) * 1.04
  xlim <- if ("xlim" %in% names(dots)) dots$xlim else c(xmin, xmax)

  ## Axis limits (y)
  ymin <- 10^-0.2
  ymax <- max(data[[varname]], wpred[[varname]], na.rm = TRUE) * 10^0.2
  ylim <- if ("ylim" %in% names(dots)) dots$ylim else c(ymin, ymax)
  data[[varname]][data[[varname]] == 0] <- ymin # set zeros to `ymin`

  ## Axis ticks (y)
  yaxis_at <- 10^(0:floor(log10(ymax)))
  yaxis_labels <- parse(text = paste0("10^", log10(yaxis_at)))

  ## Point and line styles
  points_bg_main <- "#DDDDDD"
  points_bg_ltm <- "#FFFFFF"
  points_bg_gtm <- "#882255"
  points_col_main <- "#BBBBBB"
  points_col_ltm <- "#882255"
  points_col_gtm <- "#882255"
  text_col <- "#BBBBBB"
  lines_col <- "#44AA99"
  if (inc == "cumulative") {
    data$points_bg <- points_bg_main
    data$points_col <- points_bg_main
  } else if (inc == "interval") {
    data$dt <- diff(x$time)
    m <- median(data$dt)
    dt_min <- (1 - tol) * m
    dt_max <- (1 + tol) * m
    dt_enum <- 1 + 1 * (data$dt < dt_min) + 2 * (data$dt > dt_max)
    data$points_bg <- c(points_bg_main, points_bg_ltm, points_bg_gtm)[dt_enum]
    data$points_col <- c(points_col_main, points_col_ltm, points_col_gtm)[dt_enum]
  }


  ### PLOT #############################################################

  op <- par(mar = c(5, 4, 4, 8) + 0.1, las = 1, mgp = c(3, 0.7, 0))
  on.exit(par(op))
  plot.new()
  plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i", log = "y")

  ## Axes
  box(bty = "l")
  axis(side = 1)
  axis(side = 2, at = yaxis_at, labels = yaxis_labels)

  ## Observed data
  xpd <- !any(c("xlim", "ylim") %in% names(dots))
  points(formula, data = data, xpd = xpd,
         pch = 21, bg = data$points_bg, col = data$points_col)
  if (inc == "interval") {
    text(int_inc ~ time, data = data, subset = (dt_enum != 1),
         labels = dt, pos = 3, offset = 0.3,
         cex = 0.7, font = 2, col = text_col)
  }

  ## Predicted curve
  lines(formula, data = wpred, lwd = 3, col = lines_col)

  ## Fitting window
  abline(v = x$time[c(f,l)+1], lty = 2, col = "#555555")

  ## Titles
  title(xlab = xlab)
  title(ylab = ylab)
  cstr <- x$curve
  substr(cstr, 1, 1) <- toupper(substr(cstr, 1, 1)) # capitalize first letter
  title(main = paste(cstr, "model of", inc, "incidence\n(fitted)"),
        cex.main = 0.9)

  ## Parameter estimates
  pstr1 <- paste0(names(x$theta0), " = ")
  pstr2 <- round(x$theta0, digits = 4)
  if ("K" %in% names(x$theta0)) {
    pstr2[["K"]] <- round(pstr2[["K"]])
  }
  px <- par("usr")[2] + 0.02 * diff(par("usr")[1:2])
  px <- px + max(strwidth(pstr1, cex = 0.7))
  py <- 10^(par("usr")[3] + 0.02 * diff(par("usr")[3:4]))
  text(px, py, paste(pstr1, collapse = "\n"),
       adj = c(1, 0), xpd = NA, cex = 0.7)
  text(px, py, paste(pstr2, collapse = "\n"),
       adj = c(0, 0), xpd = NA, cex = 0.7)

  ## Legend (beware: many ugly hacks here)
  lx <- par("usr")[2] + 0.02 * diff(par("usr")[1:2])
  ly <- 10^(par("usr")[4] - 0.02 * diff(par("usr")[3:4]))
  if (inc == "cumulative") {
    lstr <- c("obs", NA, NA, "pred")
    index <- c(TRUE, FALSE, FALSE, TRUE)
  } else if (inc == "interval") {
    lstr1 <- paste0("'", c("obs,", "obs,", "obs,", "pred,"), "'")
    cond <- (all(data$dt[dt_enum == 1] == m))
    rel <- c((if(cond) "=" else "~"), "<", ">", "=")
    mstr <- paste0(m, " day", if (m > 1) "s" else "")
    lstr2 <- paste0("'t ", rel, " ", mstr, "'")
    lstr <- parse(text = paste(lstr1, "~ Delta *", lstr2))
    index <- c(TRUE, any(dt_enum == 2), any(dt_enum == 3), TRUE)
  }
  legend(x = lx, y = ly,
         xpd = NA, bty = "n", cex = 0.7, seg.len = 1,
         legend = lstr[index],
         pch = c(21, 21, 21, NA)[index],
         pt.bg = c(points_bg_main, points_bg_ltm, points_bg_gtm, NA)[index],
         lty = c(NA, NA, NA, 1)[index],
         lwd = c(NA, NA, NA, 3)[index],
         col = c(points_col_main, points_col_ltm, points_col_gtm, lines_col)[index])

  invisible(NULL)
}
