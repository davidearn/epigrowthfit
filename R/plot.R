#' Plot models of epidemic growth
#'
#' @description
#' Method for plotting objects of class `"egf"`.
#'
#' @inheritParams coef.egf
#' @param type
#'   A character string, either `"cumulative"` or `"interval"`,
#'   indicating a type of incidence to plot.
#' @param xty
#'   A character string, either `"date"` or `"numeric"`,
#'   indicating how time is displayed on the bottom axis
#'   (calendar or number of days).
#' @param log
#'   A logical scalar. If `TRUE`, then incidence is measured
#'   on a logarithmic scale.
#' @param add
#'   A logical scalar. If `TRUE`, then plot elements are added
#'   to the current graphics device, which is assumed to have
#'   been started by a previous invocation of the plot method.
#' @param legend
#'   A logical scalar. If `TRUE`, then a legend is displayed
#'   in the right margin.
#' @param tol
#'   A non-negative number, used only if `inc = "interval"`.
#'   `cases[-1]` is highlighted according to `style$points_short`
#'   if `diff(date) < (1-tol)*m` and according to `style$points_long`
#'   if `diff(date) > (1+tol)*m`, where `m = median(diff(date))`.
#'   In both cases, the value of `diff(date)` is printed above
#'   the point. Assign 0 to ensure that all deviations from `m`
#'   are highlighted. Assign `Inf` to disable highlighting.
#' @param style
#'   A list of lists defining the appearance of various plot elements.
#'   See Details.
#' @param ...
#'   Optional arguments specifying additional graphical parameters.
#'   Currently, only `xlim`, `ylim`, `xlab`, `ylab`, and `main`
#'   are implemented. `xlim` can be numeric, Date, or character
#'   coercible to Date with `as.Date(xlim)`. See [graphics::plot()]
#'   and [graphics::par()] for a catalogue of graphical parameters.
#'
#' @return
#' `NULL` (invisibly).
#'
#' `get_style_default()` returns a list of lists specifying the
#' default appearance of all modifiable plot elements.
#'
#' @details
#' ## 1. Plot elements
#'
#' *Below, `date`, `cases`, `first`, `last`, `theta_init`, and
#' `theta_fit` refer to the so-named elements of `x` or `x$data`.*
#'
#' If `xty = "date"`, then the bottom axis is a calendar with ticks at
#' equally spaced Dates. If `xty = "numeric"`, then the bottom axis
#' measures the number of days since `date[1]`. Regardless of `xty`,
#' Dates `d` are assigned numeric user coordinates `d - date[1]`.
#'
#' The left axis measures interval or cumulative incidence (depending
#' on `inc`). When incidence is displayed on a logarithmic scale, zeros
#' are plotted as a positive number less than 1. They are therefore
#' distinguished from nonzero counts, which are always at least 1.
#'
#' Observed data, specified by `date` and either `cases[-1]` or
#' `cumsum(cases[-1])` (depending on `inc`), are plotted as points.
#' `cases[-1][i]` gives the number of cases observed between
#' `date[i]` and `date[i+1]`, while `cumsum(cases[-1])[i]` gives
#' the number observed between `date[1]` and `date[i+1]`. Both are
#' plotted at coordinate `date[i+1]-date[1]`.
#'
#' The fitting window is displayed as a shaded rectangle behind the
#' other plot elements. The left and right boundaries are `date[first]`
#' and `date[last]`.
#'
#' The incidence curve predicted by parameter estimates `theta_init`
#' or `theta_fit` (depending on `class(x)`) is displayed as a line
#' supported on grid points `wgrid = seq(date[first], date[last], by)`,
#' where `by = 1` for cumulative incidence and `by = median(diff(date))`
#' for interval incidence (ensuring that the interval incidence curve
#' has the correct scale; see below). The predicted curve is obtained
#' with `predict(x, time = as.integer(wgrid - date[first]))`.
#'
#' Careful interpretation of interval incidence is required if the
#' plotted time series is not equally spaced, because `cases`
#' roughly scales with `diff(date)`. That is, certain observations
#' may vary from the predicted curve not due to chance, but because
#' they give counts over fewer or more days than the typical
#' observation interval, namely `median(diff(date))`. Observations
#' for which `diff(date)` differs from `median(diff(date))` are
#' highlighted according to argument `tol` and labeled with the
#' value of `diff(date)`.
#'
#' If `annotate = TRUE` and `add = FALSE`, then a legend and the
#' parameter estimates `theta_init` or `theta_fit` (depending on
#' `class(x)`) are displayed in the right margin.
#'
#' The plot method for class "egf" displays, in addition,
#' the doubling time associated with the "egf" object and
#' the associated 95% confidence interval, obtained with
#' `confint(x, parm = "doubling_time", level = 0.95, method = "linear")`.
#'
#' ## 2. Customization
#'
#' The appearance of most of the plot elements can be controlled
#' using argument `style`. `style` must be a list containing some
#' subset of the elements below.
#'
#' \describe{
#'   \item{`box`}{
#'     A named list of arguments to [graphics::box()], affecting the
#'     appearance of the box drawn around the plot region. Currently,
#'     only `bty`, `lty`, `col`, and `lwd` are implemented.
#'   }
#'   \item{`xax`, `yax`}{
#'     Named lists of arguments to [graphics::axis()], affecting
#'     the appearance of the bottom and left axes. Currently, only
#'     `tcl`, `mgp2`, `col.axis`, `cex.axis`, and `font.axis` are
#'     implemented, where `mgp2` is the second component of the
#'     usual `mgp` argument. If `xty = "date"`, then the appearance
#'     of the minor (day or month) and major (month or year) bottom
#'     axes can be controlled separately by assigning vectors of
#'     length 2 to the elements of `xax`. For example, setting
#'     `xax = list(col.axis = c("black", "red"))` will make the
#'     minor axis black and major axis red.
#'   }
#'   \item{`xlab`, `ylab`, `main`}{
#'     Named lists of arguments to [graphics::title()], affecting
#'     the appearance of axis titles. Currently, only `line`, `adj`,
#'     `col.lab`, `cex.lab`, `font.lab`,
#'     `col.main`, `cex.main`, and `font.main` are implemented.
#'   }
#'   \item{`points_main`}{
#'     A named list of arguments to [graphics::points()],
#'     affecting the appearance of the observed data. Currently,
#'     only `pch`, `col`, `bg` and `cex` are implemented.
#'   }
#'   \item{`points_short`, `points_long`}{
#'     Alternatives to `points_main` used to highlight certain
#'     points when `inc = "interval"` (see argument `tol`).
#'   }
#'   \item{`lines`}{
#'     A named list of arguments to [graphics::lines()],
#'     affecting the appearance of the predicted incidence curve.
#'     Currently, only `lty`, `lwd`, and `col` are implemented.
#'   }
#'   \item{`window`, `confband`}{
#'     Named lists of arguments to [graphics::polygon()],
#'     affecting the appearance of the fitting window and confidence
#'     bands. Currently, only `col` and `border` are implemented.
#'   }
#'   \item{`text_hl`, `text_dbl`}{
#'     Named lists of arguments to [graphics::text()],
#'     affecting the appearance of text above highlighted points
#'     (see argument `tol`) and text giving doubling times.
#'     Currently, only `pos`, `offset`, `col`, `cex`, and `font`
#'     are implemented. `text_dbl` can further specify alignment
#'     `adj` and coordinates `x` and `y`. `x` can be numeric,
#'     Date, or character coercible to Date with `as.Date(x)`.
#'     `y` must be numeric.
#'   }
#' }
#'
#' List elements not specified by `style` are taken from
#' `get_style_default()`.
#'
#' Assigning `NULL` to any of the above possible `style` elements
#' (instead of a named list) will suppress the corresponding plot
#' element. Setting `points_main = NULL` will not suppress points
#' styled according to `points_short` and `points_long`. Hence,
#' to suppress all points, it may be necessary
#' to also set `points_short = NULL` and `points_long = NULL`.
#'
#' @seealso [egf_init()], [egf()]
#' @name plot.egf
NULL

#' @rdname plot.egf
#' @export
#' @importFrom graphics plot
plot.egf_init <- function(x, inc = "interval", xty = "date", log = TRUE,
                          add = FALSE, annotate = FALSE, tol = 0,
                          style = get_style_default(), ...) {
  ## Reuse `plot.egf()` machinery
  names(x) <- sub("theta_init", "theta_fit", names(x), fixed = TRUE)
  x <- structure(x, class = c("egf", "list"), init_flag = TRUE)
  plot(x, inc = inc, xty = xty, log = log, add = add,
       annotate = annotate, tol = tol, style = style, ...)
  invisible(NULL)
}

#' @rdname plot.egf
#' @export
#' @import graphics
#' @importFrom stats median setNames
plot.egf <- function(x, inc = "interval", xty = "date", log = TRUE,
                     add = FALSE, annotate = FALSE, tol = 0,
                     style = get_style_default(), ...) {
  ## Optional graphical parameters
  dots <- list(...)

  check(inc,
    what = "character",
    len = 1,
    opt = c("interval", "cumulative"),
    "`inc` must be one of \"interval\", \"cumulative\"."
  )
  if (inc == "cumulative") {
    check(x$data$cases[-1],
      no = anyNA,
      "Cannot calculate cumulative incidence due to missing values\n",
      "in `cases[-1]`."
    )
  }
  check(xty,
    what = "character",
    len = 1,
    opt = c("date", "numeric"),
    "`xty` must be one of \"date\", \"numeric\"."
  )
  for (a in c("log", "add", "annotate")) {
    a_val <- get(a, inherits = FALSE)
    check(a_val,
      what = "logical",
      len = 1,
      no = is.na,
      sprintf("`%s` must be TRUE or FALSE.", a)
    )
  }
  if (inc == "interval") {
    check(tol,
      what = "numeric",
      len = 1,
      val = c(0, Inf),
      no = is.na,
      "`tol` must be a non-negative number."
    )
  }
  check(style,
    what = "list",
    no = function(x) is.null(names(x)),
    "`style` must be a named list."
  )


  ### SET UP ###########################################################

  ## Flag for "egf_init" objects
  init_flag <- !is.null(attr(x, "init_flag"))

  ## Convenience
  date <- x$data$date
  cases <- x$data$cases
  i1 <- x$first
  i2 <- x$last
  d0 <- date[1]
  d1 <- date[i1]
  d2 <- date[i2]
  t1 <- days(d1, since = d0)
  t2 <- days(d2, since = d0)
  dt <- ddiff(date)

  ## Observed data
  data <- data.frame(
    time = days(date, since = d0),
    cum_inc = cumsum(c(0, cases[-1])),
    int_inc = cases,
    dt = c(NA, dt)
  )
  if (add) {
    dindex <- (data$time >= t1 - 3 & data$time <= t2 + 3)
  } else {
    dindex <- TRUE
  }

  ## A way to avoid conditional `if (inc == ...) expr1 else expr2`
  varname <- substr(inc, start = 1, stop = 3)
  varname <- paste0(varname, "_inc")
  formula <- as.formula(paste(varname, "~ time"))

  ## A way to artificially include zeros on logarithmic scale
  ymax <- max(data[[varname]], na.rm = TRUE)
  zero <- if (log) ymax^-0.04 else 0
  data[[varname]][data[[varname]] == 0] <- zero

  ## Predicted curve with confidence band
  m <- median(dt)
  wgrid <- seq(t1, t2, by = if (inc == "interval") m else 1)
  wpred <- predict(x, wgrid - t1, se = !init_flag)
  if (init_flag) {
    if (i1 > 1) {
      wpred$cum_inc <- sum(cases[2:i1]) + wpred$cum_inc
    }
    wband <- data.frame(time = wgrid, estimate = wpred[[varname]])
  } else {
    wband <- confint(wpred, level = 0.95)
    if (i1 > 1) {
      wband$cum_inc[, -1] <- sum(cases[2:i1]) + wband$cum_inc[, -1]
    }
    wband <- wband[[varname]]
    wband$time <- wgrid
  }

  ## Axis title (x)
  if (is.null(dots$xlab)) {
    if (xty == "date") {
      xlab <- "date"
    } else if (xty == "numeric") {
      xlab <- paste("days since", as.character(d0))
    }
  } else {
    xlab <- dots$xlab
  }

  ## Axis title (y)
  if (is.null(dots$ylab)) {
    ylab <- paste(inc, "incidence")
  } else {
    ylab <- dots$ylab
  }

  ## Axis title (main)
  if (is.null(dots$main)) {
    cstr <- x$curve
    substr(cstr, 1, 1) <- toupper(substr(cstr, 1, 1))
    paren <- if (init_flag) "initial guess" else "fitted"
    main <- sprintf(
      "%s model of %s incidence\n(%s)",
      cstr, inc, paren
    )
  } else {
    main <- dots$main
  }

  ## Axis limits (x)
  if (is.null(dots$xlim)) {
    xmin <- 0
    xmax <- max(data$time) * 1.04
    xlim <- c(xmin, xmax)
  } else {
    xlim <- dots$xlim
    if (is.character(xlim)) {
      xlim <- as.Date(xlim)
    }
    if (inherits(xlim, "Date")) {
      xlim <- days(xlim, since = d0)
    }
  }

  ## Axis limits (y)
  if (is.null(dots$ylim)) {
    ymin <- zero
    ymax <- if (log) ymax^1.04 else ymax * 1.04
    ylim <- c(ymin, ymax)
  } else {
    ylim <- dots$ylim
  }

  ## Style
  s <- get_style_default()
  for (pe in names(s)) {
    if (!pe %in% names(style)) {
      next
    } else if (is.null(style[[pe]])) {
      s[pe] <- list(NULL)
    } else {
      gp <- intersect(names(s[[pe]]), names(style[[pe]]))
      s[[pe]][gp] <- style[[pe]][gp]
    }
  }
  style <- s

  ## Points have sub-styles (for interval incidence,
  ## we want appearance to depend on observation interval)
  dt_min <- (1 - tol) * m
  dt_max <- (1 + tol) * m
  dt_enum <- 1 + (inc == "interval") *
    (1 * (dt < dt_min) + 2 * (dt > dt_max))
  ptys <- c("main", "short", "long")
  data$pty <- c(NA, ptys[dt_enum])


  ### PLOT #############################################################

  if (add) {
    op <- par(get("par.egf", envir = .epigrowthfit))
    on.exit(par(op))
  } else {
    op <- par(
      mar = c(4, 5, 2.7, 0.5 + 6 * annotate) + 0.1,
      bty = "l",
      xaxs = "i",
      yaxs = "i",
      las = 1
    )
    on.exit({
      assign("par.egf", par("mar", "plt"), envir = .epigrowthfit)
      assign("date.egf", date, envir = .epigrowthfit)
      par(op)
    })
    plot.new()
    plot.window(xlim = xlim, ylim = ylim, log = if (log) "y" else "")
  }

  ## Fitting window
  s <- style$window
  if (!is.null(s)) {
    l <- list(
      x = c(t1, t2, t2, t1),
      y = ylim[c(1, 1, 2, 2)]
    )
    do.call(polygon, c(l, s))
  }

  if (!add) {
    ## Box
    s <- style$box
    if (!is.null(s)) {
      do.call(box, s)
    }

    ## Axis (x)
    s <- style$xax
    if (!is.null(s)) {
      if (xty == "date") {
        l <- list(
          left = par("usr")[1],
          right = par("usr")[2],
          refdate = d0
        )
        do.call(daxis, c(l, s))
      } else if (xty == "numeric") {
        l <- list(side = 1)
        s$mgp <- c(3, s$mgp2, 0)
        s$mgp2 <- NULL
        s$
        do.call(axis, c(l, s))
      }
    }

    ## Axis (y)
    s <- style$yax
    if (!is.null(s)) {
      yax_at <- axTicks(side = 2)
      if (max(yax_at) < 1e05) {
        yax_labels <- TRUE
        long_yax_labels_flag <- FALSE
      } else {
        yax_labels <- get_labels(yax_at)
        mlw <- max(strwidth(yax_labels, units = "inches", cex = 0.85))
        long_yax_labels_flag <- (mlw / par("csi") + s$mgp2 > 3.75)
      }
      l <- list(
        side = 2,
        at = yax_at,
        labels = yax_labels
      )
      s$mgp <- c(3, s$mgp2, 0)
      s$mgp2 <- NULL
      if (is.na(s$cex.axis)) {
        s$cex.axis <- (1 - 0.25 * long_yax_labels_flag) * 0.85
      }
      do.call(axis, c(l, s))
    }
  }

  ## Observed data
  for (pty in ptys) {
    s <- style[[paste0("points_", pty)]]
    if (!is.null(s)) {
      l <- list(
        formula = formula,
        data = data,
        subset = (data$pty == pty) & dindex
      )
      do.call(points, c(l, s))
    }
  }

  ## Annotation above exceptional points
  s <- style$text_hl
  if (!is.null(s)) {
    l <- list(
      formula = int_inc ~ time,
      data = data,
      labels = data$dt,
      subset = c(NA, (data$pty != ptys[1] & dindex)[-1])
    )
    do.call(text, c(l, s))
  }

  ## Confidence band
  if (!init_flag) {
    s <- style$confband
    if (!is.null(s)) {
      l <- list(
        x = c(wband$time, rev(wband$time)),
        y = c(wband$lower, rev(wband$upper))
      )
      do.call(polygon, c(l, s))
    }
  }

  ## Predicted curve
  s <- style$lines
  if (!is.null(s)) {
    l <- list(formula = estimate ~ time, data = wband)
    do.call(lines, c(l, s))
  }

  ## Doubling time
  s <- style$text_dbl
  if (!is.null(s) && !init_flag && inc == "interval") {
    ci <- confint(x, parm = "doubling_time", level = 0.95, trace = FALSE)
    dblstr <- sprintf(
      "doubling time:\n%.1f (%.1f, %.1f) days",
      ci[1], ci[2], ci[3]
    )
    wrange <- range(wband$estimate, na.rm = TRUE)
    if (is.na(s$y)) {
      if (log) {
        s$y <- wrange[1] * 10^(0.25 * diff(log10(wrange)))
      } else {
        s$y <- wrange[1] + 0.25 * diff(wrange)
      }
    }
    if (is.character(s$x)) {
      s$x <- as.Date(s$x)
    }
    if (inherits(s$x, "Date")) {
      s$x <- days(s$x, since = d0)
    }
    if (is.na(s$x)) {
      s$x <- min(wband$time[wband$estimate > s$y], na.rm = TRUE)
    }
    l <- list(labels = dblstr, xpd = NA)
    do.call(text, c(l, s))
  }

  if (!add) {
    ## Axis title (x)
    s <- style$xlab
    if (!is.null(s)) {
      l <- list(xlab = xlab)
      do.call(title, c(l, s))
    }

    ## Axis title (y)
    s <- style$ylab
    if (!is.null(s)) {
      l <- list(ylab = ylab)
      do.call(title, c(l, s))
    }

    ## Axis title (main)
    s <- style$main
    if (!is.null(s)) {
      l <- list(main = main)
      do.call(title, c(l, s))
    }

    if (annotate) {
      ## Parameter estimates
      pstr <- sprintf("%s = ", names(x$theta_fit))
      pvec <- round(x$theta_fit,
        digits = ifelse(names(x$theta_fit) == "K", 0, 4)
      )
      px <- par("usr")[2] + 0.02 * diff(par("usr")[1:2]) +
        max(strwidth(pstr, cex = 0.7))
      py <- par("usr")[3] + 0.02 * diff(par("usr")[3:4])
      if (log) {
        py <- 10^py
      }
      text(px, py, paste(pstr, collapse = "\n"),
           adj = c(1, 0), xpd = NA, cex = 0.7)
      text(px, py, paste(pvec, collapse = "\n"),
           adj = c(0, 0), xpd = NA, cex = 0.7)

      ## Legend
      lx <- par("usr")[2] + 0.02 * diff(par("usr")[1:2])
      ly <- par("usr")[4] - 0.02 * diff(par("usr")[3:4])
      if (log) {
        ly <- 10^ly
      }
      if (inc == "cumulative") {
        lstr <- c("obs", NA, NA, "pred")
        index <- c(TRUE, FALSE, FALSE, TRUE)
      } else if (inc == "interval") {
        cond <- all(dt[dt_enum == 1] == m)
        lstr <- sprintf(
          "'%s,' ~ Delta * 't %s %g day%s'",
          c("obs", "obs", "obs", "pred"),
          c(if (cond) "=" else "~", "<", ">", "="),
          m,
          if (m > 1) "s" else ""
        )
        lstr <- parse(text = lstr)
        index <- c(TRUE, any(dt_enum == 2), any(dt_enum == 3), TRUE)
      }
      legend(x = lx, y = ly,
        legend = lstr[index],
        xpd = NA, bty = "n", cex = 0.7, seg.len = 1,
        pch = c(style$points_main$pch, style$points_short$pch, style$points_long$pch, NA)[index],
        pt.bg = c(style$points_main$bg, style$points_short$bg, style$points_long$bg, NA)[index],
        lty = c(NA, NA, NA, style$lines$lty)[index],
        lwd = c(NA, NA, NA, style$lines$lwd)[index],
        col = c(style$points_main$col, style$points_short$col, style$points_long$col, style$lines$col)[index]
      )
    }
  }

  invisible(NULL)
}

#' @rdname plot.egf
#' @export
get_style_default <- function() {
  list(
    box = list(bty = "l", lty = 1, lwd = 1, col = "black"),
    xax = list(tcl = -0.2, mgp2 = c(0.25, 1.25), col.axis = "black", cex.axis = c(0.85, 1.15), font.axis = 1),
    yax = list(tcl = -0.5, mgp2 = 0.7, col.axis = "black", cex.axis = NA, font.axis = 1),
    xlab = list(line = 3, adj = 0.5, col.lab = "black", cex.lab = 1, font.lab = 1),
    ylab = list(line = 4, adj = 0.5, col.lab = "black", cex.lab = 1, font.lab = 1),
    main = list(line = 1, adj = 0.5, col.main = "black", cex.main = 1, font.main = 2),
    points_main = list(pch = 21, col = "#BBBBBB", bg = "#DDDDDD", cex = 1),
    points_short = list(pch = 1, col = "#882255", bg = NA, cex = 1),
    points_long = list(pch = 16, col = "#882255", bg = NA, cex = 1),
    lines = list(lty = 1, lwd = 2.5, col = "#44AA99"),
    window = list(col = "#DDCC7740", border = NA),
    confband = list(col = "#44AA9960", border = NA),
    text_hl = list(pos = 3, offset = 0.3, col = "#BBBBBB", cex = 0.7, font = 2),
    text_dbl = list(x = NA, y = NA, adj = c(0, 0.5), pos = 4, offset = 1, col = "black", cex = 0.7, font = 1)
  )
}

