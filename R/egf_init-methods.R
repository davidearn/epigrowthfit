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
#' @param add A logical scalar. If `TRUE`, then the fitting window
#'   and predicted curve (and nothing else) are added to the current
#'   graphics device.
#' @param annotate A logical scalar. If `TRUE`, then a legend and
#'   a list of parameter estimates are added to the right margin.
#'   Ignored if `add = TRUE`.
#' @param tol A non-negative number used only if `inc = "interval"`.
#'   `cases[i]` is highlighted according to `point_style_short` if
#'   `diff(time)[i] < (1-tol)*m` and according to `point_style_long`
#'   if `diff(time)[i] > (1+tol)*m`, where `m = median(diff(time))`.
#'   In both cases, the value of `diff(time)[i]` is printed above
#'   the point. Assign 0 to ensure that all deviations from `m` are
#'   highlighted. Assign `Inf` to disable highlighting.
#' @param polygon_style A named list of arguments to
#'   [`polygon()`][graphics::polygon()], affecting the appearance
#'   of the fitting window. Currently, only `col` and `border` are
#'   implemented.
#' @param line_style A named list of arguments to
#'   [`lines()`][graphics::lines()], affecting the appearance
#'   of the predicted curve. Currently, only `lty`, `lwd`, and `col`
#'   are implemented.
#' @param point_style_main A named list of arguments to
#'   [`points()`][graphics::points()], affecting the appearance
#'   of the observed data. Currently, only `pch`, `col`, `bg` and `cex`
#'   are implemented.
#' @param point_style_short,point_style_long
#'   Alternatives to `point_style_main` used to highlight
#'   certain points when `inc = "interval"`. See argument `tol`.
#' @param text_style A named list of arguments to
#'   [`text()`][graphics::text()], affecting the appearance
#'   of text printed above highlighted points. See argument `tol`.
#'   Currently, only `pos`, `offset`, `col`, `cex`, and `font`
#'   are implemented.
#' @param ... Optional arguments. Used only by the `plot` method
#'   to specify graphical parameters. See [`plot()`][base::plot()]
#'   and [`par()`][graphics::par()]. Currently, only `xlim`, `ylim`,
#'   `xlab`, `ylab`, and `main` are implemented. Any additional
#'   parameters will be ignored.
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
#'     Equal to `object$eval_cum_inc(time)`.
#'   }
#'   \item{`int_inc`}{Expected interval incidence given interval
#'     endpoints `time`, conditional on parameter vector
#'     `object$theta0`. Equal to `diff(object$eval_cum_inc(time))`
#'     if `length(time) >= 2` and omitted otherwise.
#'   }
#' }
#'
#' The `plot` method returns `NULL` invisibly.
#'
#' @details
#' ## Plot elements
#'
#' The bottom axis measures the number of days since `x$date[1]`. The
#' left axis measures interval or cumulative incidence (depending on
#' `inc`) on a log scale. Zeros are plotted as if they were `10^-0.2`,
#' and are therefore distinguished from nonzero counts, which are
#' always at least 1.
#'
#' Observed data, specified by `x$time` and either `x$cases` or
#' `cumsum(x$cases)` (depending on `inc`), are plotted as points.
#' `cases[i]` gives the number of cases observed between `time[i]`
#' and `time[i+1]`, and `cumsum(cases)[i]` the number observed
#' between `time[1]` and `time[i+1]`. Both are plotted at `time[i+1]`.
#'
#' The fitting window, specified by indices `x$first` and `x$last`, is
#' displayed as a shaded rectangle behind the other plot elements. The
#' left and right boundaries occur at `time[first]` and `time[last+1]`.
#'
#' The incidence curve predicted by initial parameter estimates
#' `x$theta0` is displayed as a line supported on grid points
#' `wgrid = seq(time[first], time[last+1], by)`, where `by = 1`
#' for cumulative incidence and `by = median(diff(time))` for
#' interval incidence (to ensure that the curve has the correct
#' scale; see below). The predicted curve is obtained with
#' `predict(x, wgrid)`.
#'
#' Careful interpretation of interval incidence is required if
#' the plotted time series is not equally spaced, because `cases`
#' roughly scales with `diff(time)`. That is, certain observations
#' may vary from the predicted curve not due to chance, but because
#' they represent a count over fewer or more days than the typical
#' observation interval, namely `median(diff(time))`. Observations
#' for which `diff(time)` differs from `median(diff(time))` are
#' highlighted according to argument `tol` and labeled with the
#' value of `diff(time)`.
#'
#' If `add = FALSE` and `annotate = TRUE`, then a legend and the
#' initial parameter estimates `x$theta0` are printed in the right
#' margin.
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
  first <- x$first
  last <- x$last
  cat("Pass this \"egf_init\" object to `egf()` to fit", cstr, "\n")
  cat(bstr, dstr, "\n")
  cat("\n")
  cat("Fitting window:\n")
  cat("\n")
  cat("index   ", first, ":", last, "\n", sep = "")
  cat(" date   (", as.character(x$date[first]), ", ", as.character(x$date[last+1]), "]\n", sep = "")
  cat("cases   ", sum(x$cases[first:last]), " of ", sum(x$cases), "\n", sep = "")
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
    cum_inc = object$eval_cum_inc(time)
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
#' @importFrom scales alpha
plot.egf_init <- function(x, inc = "interval",
                          add = FALSE, annotate = TRUE, tol = 0,
                          polygon_style = list(col = "#DDCC7740", border = NA),
                          line_style = list(lty = 1, lwd = 3, col = "#44AA99"),
                          point_style_main = list(pch = 21, col = "#BBBBBB", bg = "#DDDDDD", cex = 1),
                          point_style_short = list(pch = 1, col = "#882255", bg = NA, cex = 1),
                          point_style_long = list(pch = 16, col = "#882255", bg = NA, cex = 1),
                          text_style = list(pos = 3, offset = 0.3, col = "#BBBBBB", cex = 0.7, font = 2),
                          ...) {
  if (!is.character(inc) || length(inc) != 1 ||
        !inc %in% c("interval", "cumulative")) {
    stop("`inc` must be one of \"interval\", \"cumulative\".")
  }
  if (!is.logical(add) || length(add) != 1 || is.na(add)) {
    stop("`add` must be `TRUE` or `FALSE`.")
  }
  if (!is.logical(annotate) || length(annotate) != 1 || is.na(annotate)) {
    stop("`annotate` must be `TRUE` or `FALSE`.")
  }
  if (inc == "interval") {
    if (!is.numeric(tol) || length(tol) != 1 || !isTRUE(tol >= 0)) {
      stop("`tol` must be a non-negative number.")
    }
  }
  l <- list(polygon_style, line_style, text_style,
            point_style_main, point_style_short, point_style_long)
  if (!all(sapply(l, is.list))) {
    stop("All \"_style\" arguments must be lists.")
  }


  ### SET UP ###########################################################

  ## Optional graphical parameters
  dots <- list(...)

  ## Observed data
  data <- data.frame(
    time = x$time[-1],
    cum_inc = cumsum(x$cases),
    int_inc = x$cases,
    dt = diff(x$time)
  )
  m <- median(data$dt)

  ## Predicted curve
  tf <- x$time[x$first]
  tl <- x$time[x$last+1]
  wgrid <- seq(tf, tl, by = if (inc == "interval") m else 1)
  wpred <- predict(x, wgrid)[c("time", "cum_inc", "int_inc")]
  wpred$int_inc <- c(NA, wpred$int_inc)

  ## A way to avoid conditional `if (inc = ...) ... else ...`
  varname <- substr(inc, start = 1, stop = 3) # first three characters
  varname <- paste0(varname, "_inc")
  formula <- as.formula(paste(varname, "~ time"))

  ## Titles
  xlab <- if ("xlab" %in% names(dots)) dots$xlab else paste("days since", as.character(x$date[1]))
  ylab <- if ("ylab" %in% names(dots)) dots$ylab else paste(inc, "incidence")
  cstr <- x$curve
  substr(cstr, 1, 1) <- toupper(substr(cstr, 1, 1)) # capitalize first letter
  main <- if ("main" %in% names(dots)) dots$main else paste(cstr, "model of", inc, "incidence\n(initial guess)")

  ## Axis limits (x)
  xmin <- 0
  xmax <- max(x$time) * 1.04
  xlim <- if ("xlim" %in% names(dots)) dots$xlim else c(xmin, xmax)

  ## Axis limits (y)
  ymin <- 10^-0.2
  ymax <- max(data[[varname]]) * 10^0.2
  ylim <- if ("ylim" %in% names(dots)) dots$ylim else c(ymin, ymax)
  data[[varname]][data[[varname]] == 0] <- ymin # set zeros to `ymin`

  ## Axis ticks (y)
  yaxis_at <- 10^(0:floor(log10(ymax)))
  yaxis_labels <- parse(text = paste0("10^", log10(yaxis_at)))

  ## Styles
  for (a in grep("_style", names(formals(plot.egf_init)), value = TRUE)) {
    l1 <- eval(formals(plot.egf_init)[[a]]) # default style
    l2 <- get(a) # passed style
    inter <- intersect(names(l1), names(l2))
    l1[inter] <- l2[inter]
    assign(a, l1)
  }

  ## Style for each point
  ## (for interval incidence, style depends on observation interval)
  dt_min <- (1 - tol) * m
  dt_max <- (1 + tol) * m
  dt_enum <- 1 + (inc == "interval") * (1 * (data$dt < dt_min) + 2 * (data$dt > dt_max))
  point_style_list <- c("main", "short", "long")
  data$point_style <- point_style_list[dt_enum]


  ### PLOT #############################################################

  if (add) {
    sp <- get("par", envir = .egf_env)
    sp$yaxp <- NULL # keeping `yaxp` causes `par()` to throw an error
    op <- par(sp)
  } else {
    op <- par(mar = c(4, 4, 2, 7 * annotate) + 0.5, mgp = c(3, 0.7, 0), las = 1)
    plot.new()
    plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i", log = "y")
  }

  ## Fitting window
  l <- list(
    x = c(tf, tl, tl, tf),
    y = c(ymin, ymin, ymax, ymax)
  )
  do.call(polygon, c(l, polygon_style))
  windex <- (data$time >= tf - 6 & data$time <= tl + 6)

  ## Axes
  if (!add) {
    box(bty = "l")
    axis(side = 1)
    axis(side = 2, at = yaxis_at, labels = yaxis_labels)
  }

  ## Observed data
  for (ps in point_style_list) {
    l <- list(
      formula = formula,
      data = data,
      subset = (data$point_style == ps) & (if (add) windex else TRUE),
      xpd = !any(c("xlim", "ylim") %in% names(dots))
    )
    do.call(points, c(l, get(paste0("point_style_", ps))))
  }

  ## Annotation above exceptional points
  if (any(dt_enum != 1)) {
    l <- list(
      formula = int_inc ~ time,
      data = data,
      subset = (dt_enum != 1) & (if (add) windex else TRUE)
    )
    do.call(text, c(l, text_style))
  }

  ## Predicted curve
  l <- list(
    formula = formula,
    data = wpred
  )
  do.call(lines, c(l, line_style))

  if (!add) {
    ## Titles
    title(xlab = xlab)
    title(ylab = ylab)
    title(main = main, cex.main = 0.9)

    if (annotate) {
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

      ## Legend (beware: some ugly hacks here)
      lx <- par("usr")[2] + 0.02 * diff(par("usr")[1:2])
      ly <- 10^(par("usr")[4] - 0.02 * diff(par("usr")[3:4]))
      if (inc == "cumulative") {
        lstr <- c("obs", NA, NA, "pred")
        index <- c(TRUE, FALSE, FALSE, TRUE)
      } else if (inc == "interval") {
        lstr1 <- paste0("'", c("obs,", "obs,", "obs,", "pred,"), "'")
        cond <- (all(data$dt[dt_enum == 1] == m))
        rel <- c(if(cond) "=" else "~", "<", ">", "=")
        mstr <- paste0(m, " day", if (m > 1) "s" else "")
        lstr2 <- paste0("'t ", rel, " ", mstr, "'")
        lstr <- parse(text = paste(lstr1, "~ Delta *", lstr2))
        index <- c(TRUE, any(dt_enum == 2), any(dt_enum == 3), TRUE)
      }
      legend(x = lx, y = ly,
             xpd = NA, bty = "n", cex = 0.7, seg.len = 1,
             legend = lstr[index],
             pch = c(point_style_main$pch, point_style_short$pch, point_style_long$pch, NA)[index],
             pt.bg = c(point_style_main$bg, point_style_short$bg, point_style_long$bg, NA)[index],
             lty = c(NA, NA, NA, line_style$lty)[index],
             lwd = c(NA, NA, NA, line_style$lwd)[index],
             col = c(point_style_main$col, point_style_short$col, point_style_long$col, line_style$col)[index])
    }
  }

  if (!add) {
    assign("par", par(no.readonly = TRUE), envir = .egf_env)
  }
  par(op)
  invisible(NULL)
}
