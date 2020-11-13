#' Plot models of epidemic growth
#'
#' @description
#' Methods for plotting objects of class "egf_init" or "egf".
#'
#' @param x
#'   An "egf_init" or "egf" object.
#' @param inc
#'   One of `"cumulative"` and `"interval"`,
#'   indicating a type of incidence to plot.
#' @param xty
#'   One of `"date"` and `"numeric"`,
#'   indicating how time is displayed on the bottom axis.
#' @param log
#'   A logical scalar. If `TRUE`, then incidence is displayed
#'   on a logarithmic scale on the left axis.
#' @param add
#'   A logical scalar. If `TRUE`, then plot elements are added
#'   to the current graphics device, which is assumed to have
#'   been started by a previous invocation of the plot method.
#' @param annotate
#'   A logical scalar. If `TRUE`, then a legend and a list of
#'   parameter estimates are displayed in the right margin.
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
#'   See Details 2.
#' @param ...
#'   Optional arguments specifying additional graphical parameters.
#'   Currently, only `xlim`, `ylim`, `xlab`, `ylab`, and `main`
#'   are implemented. `xlim` can be numeric, Date, or character
#'   coercible to Date with `as.Date(xlim)`. See [graphics::plot()]
#'   and [graphics::par()] for a catalogue of graphical parameters.
#'
#' @return
#' `plot.egf_init()` and `plot.egf()` return `NULL` (invisibly).
#'
#' `get_style_default()` returns a list of lists specifying the
#' default appearance of all (modifiable) plot elements.
#'
#' @details
#' ## 1. Plot elements
#'
#' *Below, `date`, `time`, `cases`, `first`, `last`, `theta_init`,
#' and `theta_fit` refer to the so-named elements of `x` or `x$init`.*
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
#' `cases[-1][i]` gives the number of cases observed between `date[i]`
#' and `date[i+1]`, while `cumsum(cases[-1])[i]` gives the number
#' observed between `date[1]` and `date[i+1]`. Both are plotted
#' at coordinate `date[i+1]-date[1]`.
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
#' plotted time series is not equally spaced, because `cases` roughly
#' scales with `diff(date)`. That is, certain observations may vary
#' from the predicted curve not due to chance, but because they
#' represent a count over fewer or more days than the typical
#' observation interval, namely `median(diff(date))`. Observations for
#' which `diff(date)` differs from `median(diff(date))` are highlighted
#' according to argument `tol` and labeled with the value of `diff(date)`.
#'
#' If `annotate = TRUE` and `add = FALSE`, then a legend and the
#' parameter estimates `theta_init` or `theta_fit` (depending on
#' the `class(x)`) are displayed in the right margin.
#'
#' The plot method for class "egf" displays, in addition,
#' the doubling time associated with the "egf" object and
#' the associated 95% confidence interval, obtained with
#' `compute_doubling_time(rev(unname(confint(x))))`.
#'
#' ## 2. Customization
#'
#' The appearance of most of the plot elements can be controlled
#' using argument `style`. `style` must be a list containing some
#' subset of the elements below.
#'
#' \describe{
#'   \item{`date`}{
#'     A named list of arguments to internal function `daxis()`
#'     (a subset of `tcl`, `mgp2`, `col.axis`, and `cex.axis`),
#'     affecting the appearance of the bottom axis if `xty = "date"`.
#'   }
#'   \item{`points_main`}{
#'     A named list of arguments to [graphics::points()],
#'     affecting the appearance of the observed data. Currently,
#'     only `pch`, `col`, `bg` and `cex` are implemented.
#'     Set to `NULL` to suppress.
#'   }
#'   \item{`points_short`, `points_long`}{
#'     Alternatives to `points_main` used to highlight certain
#'     points when `inc = "interval"` (see argument `tol`).
#'     Set to `NULL` to suppress.
#'
#'   }
#'   \item{`lines`}{
#'     A named list of arguments to [graphics::lines()],
#'     affecting the appearance of the predicted incidence curve.
#'     Currently, only `lty`, `lwd`, and `col` are implemented.
#'     Set to `NULL` to suppress.
#'   }
#'   \item{`window`, `confband`}{
#'     Named lists of arguments to [graphics::polygon()],
#'     affecting the appearance of the fitting window and confidence
#'     bands. Currently, only `col` and `border` are implemented.
#'     Set to `NULL` to suppress.
#'   }
#'   \item{`text_hl`, `text_dbl`}{
#'     Named lists of arguments to [graphics::text()],
#'     affecting the appearance of text above highlighted points
#'     (see argument `tol`) and text giving doubling times.
#'     Currently, only `pos`, `offset`, `col`, `cex`, and `font`
#'     are implemented. `text_dbl` can further specify coordinates
#'     `x` and `y` and alignment `adj`, and `x` can be numeric,
#'     Date, or character coercible to Date with `as.Date(x)`.
#'     Set to `NULL` to suppress.
#'   }
#' }
#'
#' List elements not specified by `style` are taken from
#' `get_style_default()`.
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
  ## Disguise `x` as an "egf" object to reuse `plot.egf()` machinery
  x <- list(
    init = x,
    theta_fit = x$theta_init,
    eval_cum_inc = x$eval_cum_inc
  )
  x <- structure(x, class = c("egf", "list"), init_flag = TRUE)
  plot(x, inc = inc, xty = xty, log = log,
       add = add, annotate = annotate, tol = tol,
       style = style, ...)
  invisible(NULL)
}

#' @rdname plot.egf
#' @export
#' @import graphics
#' @importFrom stats median
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
    check(x$init$cases[-1],
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
  date <- x$init$date
  cases <- x$init$cases
  i1 <- x$init$first
  i2 <- x$init$last
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

  ## Predicted curve
  m <- median(dt)
  wgrid <- seq(t1, t2, by = if (inc == "interval") m else 1)
  wpred <- predict(x, wgrid - t1)[c("time", "cum_inc", "int_inc")]
  wpred$time <- wgrid
  if (i1 > 1) {
    wpred$cum_inc <- sum(cases[2:i1]) + wpred$cum_inc
  }

  ## A way to avoid conditional `if (inc == ...) expr1 else expr2`
  varname <- substr(inc, start = 1, stop = 3)
  varname <- paste0(varname, "_inc")
  formula <- as.formula(paste(varname, "~ time"))

  ## A way to artificially include zeros on logarithmic scale
  ymax <- max(data[[varname]], na.rm = TRUE)
  zero <- if (log) ymax^-0.04 else 0
  data[[varname]][data[[varname]] == 0] <- zero

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
    cstr <- x$init$curve
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
      las = 1,
      mgp = c(3, 0.7, 0),
      cex.axis = 0.85,
      cex.main = 0.9
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
    box()

    ## Axis (x)
    if (xty == "date") {
      l <- list(
        left = par("usr")[1],
        right = par("usr")[2],
        refdate = d0
      )
      do.call(daxis, c(l, style$date))
    } else if (xty == "numeric") {
      axis(side = 1)
    }

    ## Axis (y)
    yax_at <- axTicks(side = 2)
    if (max(yax_at) < 1e05) {
      yax_labels <- TRUE
      long_yax_labels_flag <- FALSE
    } else {
      yax_labels <- get_labels(yax_at)
      mlw <- max(strwidth(yax_labels, units = "inches", cex = par("cex.axis")))
      long_yax_labels_flag <- (mlw / par("csi") + par("mgp")[2] > 3.75)
    }
    axis(side = 2, at = yax_at, labels = yax_labels,
         cex.axis = (1 - 0.25 * long_yax_labels_flag) * par("cex.axis"))
  }

  ## Observed data
  for (pty in ptys) {
    s <- style[[paste0("points_", pty)]]
    if (!is.null(s)) {
      l <- list(
        formula = formula,
        data = data,
        subset = (data$pty == pty) & dindex,
        xpd = is.null(dots$xlim) && is.null(dots$ylim)
      )
      do.call(points, c(l, s))
    }
  }

  ## Annotation above exceptional points
  s <- style$text_hl
  is_exceptional <- (data$pty != ptys[1] & dindex)
  if (!is.null(s) && isTRUE(any(is_exceptional))) {
    l <- list(
      formula = int_inc ~ time,
      data = data,
      labels = data$dt,
      subset = is_exceptional
    )
    do.call(text, c(l, s))
  }

  ## Predicted curve
  s <- style$lines
  if (!is.null(s)) {
    l <- list(formula = formula, data = wpred)
    do.call(lines, c(l, s))
  }

  ## Doubling time
  s <- style$text_dbl
  if (!is.null(s) && !init_flag && inc == "interval") {
    estimate <- compute_doubling_time(x)
    sink(nullfile())
    ci <- confint(x, parm = "r", level = 0.95, method = "linear")
    sink(NULL)
    ci <- rev(unname(compute_doubling_time(ci)))
    dblstr <- sprintf(
      "doubling time:\n%.1f (%.1f, %.1f) days",
      estimate, ci[1], ci[2]
    )
    wrange <- range(wpred$int_inc, na.rm = TRUE)
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
      s$x <- min(wpred$time[wpred$int_inc > s$y], na.rm = TRUE)
    }
    l <- list(labels = dblstr, xpd = NA)
    do.call(text, c(l, s))
  }

  if (!add) {
    ## Titles
    title(xlab = xlab, line = 3)
    title(ylab = ylab, line = 4)
    title(main = main, line = 1)

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
#' @keywords internal
#' @export
get_style_default <- function() {
  list(
    date = list(tcl = -0.2, mgp2 = c(0.05, 1), col.axis = c("black", "black"), cex.axis = c(0.7, 0.85)),
    points_main = list(pch = 21, col = "#BBBBBB", bg = "#DDDDDD", cex = 1),
    points_short = list(pch = 1, col = "#882255", bg = NA, cex = 1),
    points_long = list(pch = 16, col = "#882255", bg = NA, cex = 1),
    lines = list(lty = 1, lwd = 3, col = "#44AA99"),
    window = list(col = "#DDCC7740", border = NA),
    confband = list(col = "#44AA9940", border = NA),
    text_hl = list(pos = 3, offset = 0.3, col = "#BBBBBB", cex = 0.7, font = 2),
    text_dbl = list(x = NA, y = NA, adj = c(0, 0.5), pos = NULL, offset = 1, col = "black", cex = 0.7, font = 1)
  )
}

#' Plot simulations
#'
#' @description
#' A method for plotting simulated incidence curves specified
#' by objects of class "egf_sim".
#'
#' @param x
#'   An "egf_init" or "egf" object.
#' @param inc
#'   One of `"cumulative"` and `"interval"`,
#'   indicating a type of incidence to plot.
#' @param col_sim,col_pred
#'   Character or numeric scalars specifying colours
#'   for simulated and predicted incidence curves, respectively.
#' @param ...
#'   Optional arguments specifying additional graphical parameters.
#'   Currently, only `xlim`, `ylim`, `xlab`, `ylab`, and `main` are
#'   implemented. See [graphics::plot()] and [graphics::par()] for
#'   a catalogue of graphical parameters.
#'
#' @return
#' `NULL` (invisibly).
#'
#' @details
#' ## Plot elements
#'
#' The bottom axis measures the number of days since
#' `with(x$object$init, date[first])`.
#' The left axis measures interval or cumulative incidence
#' (depending on `inc`).
#'
#' The predicted incidence curve is obtained as `pred$cum_inc` or
#' `pred$int_inc`, where `pred = predict(x$object, time = x$time)`.
#'
#' The simulated incidence curves are obtained as the columns of
#' `x$cum_inc` or `x$int_inc`. These are plotted together behind
#' the predicted incidence curve.
#'
#' To help visualize the distribution of simulated incidence at a
#' given time, assign `col_sim` a sufficiently transparent colour.
#'
#' @seealso [simulate.egf()]
#' @export
#' @import graphics
plot.egf_sim <- function(x, inc = "cumulative",
                         col_pred = "#44AA99", col_sim = "#BBBBBB66", ...) {
  check(inc,
    what = "character",
    len = 1,
    opt = c("interval", "cumulative"),
    "`inc` must be one of \"interval\", \"cumulative\"."
  )


  ### SET UP ###########################################################

  ## Optional graphical parameters
  dots <- list(...)

  ## Predicted curve
  pred <- predict(x$object, time = x$time)[c("time", "cum_inc", "int_inc")]

  ## A way to avoid conditional `if (inc = ...) expr1 else expr2`
  varname <- substr(inc, start = 1, stop = 3)
  varname <- paste0(varname, "_inc")

  ## Number of simulations
  nsim <- ncol(x[[varname]])

  ## Axis title (x)
  if (is.null(dots$xlab)) {
    xlab <- paste("days since", as.character(x$refdate))
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
    cstr <- x$object$init$curve
    substr(cstr, 1, 1) <- toupper(substr(cstr, 1, 1))
    main <- sprintf(
      "%s model of %s incidence\n(%d simulations)",
      cstr, inc, nsim
    )
  } else {
    main <- dots$main
  }

  ## Axis limits (x)
  if (is.null(dots$xlim)) {
    xmin <- min(x$time)
    xmax <- max(x$time) * 1.04
    xlim <- c(xmin, xmax)
  } else {
    xlim <- dots$xlim
  }

  ## Axis limits (y)
  if (is.null(dots$ylim)) {
    ymin <- 0
    ymax <- max(x[[varname]], na.rm = TRUE) * 1.04
    ylim <- c(ymin, ymax)
  } else {
    ylim <- dots$ylim
  }


  ### PLOT #############################################################

  op <- par(
    mar = c(4, 5, 2.7, 0.5) + 0.1,
    bty = "l",
    xaxs = "i",
    yaxs = "i",
    las = 1,
    mgp = c(3, 0.7, 0),
    cex.axis = 0.85,
    cex.main = 0.9
  )
  on.exit(par(op))
  plot.new()
  plot.window(xlim = xlim, ylim = ylim)

  ## Simulated data
  for (j in 1:nsim) {
    lines(x$time, x[[varname]][, j], lwd = 2, col = col_sim)
  }

  ## Predicted curves
  lines(pred$time, pred[[varname]], lwd = 3, col = col_pred)

  ## Box
  box()

  ## Axis (x)
  axis(side = 1)

  ## Axis (y)
  yax_at <- axTicks(side = 2)
  if (max(yax_at) < 1e05) {
    yax_labels <- TRUE
    long_yax_labels_flag <- FALSE
  } else {
    yax_labels <- get_labels(yax_at)
    mlw <- max(strwidth(yax_labels, units = "inches", cex = par("cex.axis")))
    long_yax_labels_flag <- (mlw / par("csi") + par("mgp")[2] > 3.75)
  }
  axis(side = 2, at = yax_at, labels = yax_labels,
       cex.axis = (1 - 0.25 * long_yax_labels_flag) * par("cex.axis"))

  ## Titles
  title(xlab = xlab, line = 3)
  title(ylab = ylab, line = 4)
  title(main = main, line = 1)

  invisible(NULL)
}

#' Plot smoothed incidence data
#'
#' @description
#' A method for plotting objects of class "smooth_cases".
#'
#' @param x
#'   A "smooth_cases" object.
#' @param v
#'   An integer vector indexing `x$date`, indicating that
#'   vertical lines should be drawn at dates `x$date[v]`.
#'   Alternatively, a Date vector with elements between
#'   `min(x$date)` and `max(x$date)` or a character vector
#'   coercible to such a Date vector via `as.Date()`
#'   (e.g., "YYYY-MM-DD"). In this case, coercion from Date
#'   to index is done by `which.min(abs(x$date - v))`.
#' @param ...
#'   For `plot.smooth_cases()`, unused optional arguments.
#'   For `dline()`, optional arguments to [graphics::abline()],
#'   such as `lty`, `lwd`, and `col`. These are recycled
#'   to the length of `v`.
#'
#' @return
#' `NULL` (invisibly).
#'
#' @details
#' ## Plot elements
#'
#' The bottom axis is a calendar with ticks at equally spaced dates.
#' Dates `d` have numeric user coordinates `d - x$date[1]`. The left
#' axis measures interval incidence, log(1+x)-transformed if
#' `x$log = TRUE`.
#'
#' Plotted as points are the interval incidence data specified by
#' `x$date` and `x$cases`. Plotted as a solid line is the cubic
#' spline fit to the (possibly transformed) data specified by `x$ss`.
#'
#' Dotted lines are drawn at peaks (red) and troughs (blue) in
#' the spline. They are labeled in the top margin with the index
#' of `x$date` corresponding to the approximate date of the peak
#' or trough in the spline.
#'
#' @export
#' @import graphics
plot.smooth_cases <- function(x, ...) {
  ### SET UP ###########################################################

  ## Observed data
  data <- data.frame(
    time = days(x$date, since = x$date[1]),
    cases = x$cases,
    log_cases = log10(1 + x$cases)
  )
  varname <- if (x$log) "log_cases" else "cases"
  formula <- as.formula(paste(varname, "~ time"))

  ## Times of peaks and troughs
  time_peaks <- data$time[x$peaks]
  time_troughs <- data$time[x$troughs]

  ## Axis limits
  xlim <- range(data$time)
  ymax <- max(data[[varname]], na.rm = TRUE) * 1.04
  ylim <- c(0, ymax)

  ## Axis titles
  xlab <- "date"
  ylab <- if (x$log) expression(log[10] * "(1+cases)") else "cases"

  ## Colour palette
  col_points <- "#CCCCCC"
  col_peaks <- "#BB5566"
  col_troughs <- "#004488"

  ## Lines for tick labels in top margin
  mgp2_peaks <- if (length(x$troughs) > 0) 0.6 else 0
  mgp2_troughs <- 0


  ### PLOT #############################################################

  op <- par(
    mar = c(2.9, 4.1, 3.1, 0.1) + 0.5,
    bty = "l",
    xaxs = "i",
    yaxs = "i"
  )
  on.exit({
    assign("par.smooth_cases", par("mar", "plt"), envir = .epigrowthfit)
    assign("date.smooth_cases", x$date, envir = .epigrowthfit)
    par(op)
  })

  plot.new()
  plot.window(xlim = xlim, ylim = ylim)
  points(formula, data = data, col = col_points)
  lines(y ~ x, data = predict(x$ss, data$time), lwd = 2)
  abline(v = time_peaks, lty = 3, col = col_peaks)
  abline(v = time_troughs, lty = 3, col = col_troughs)
  box()
  daxis(left = par("usr")[1], right = par("usr")[2],
        refdate = x$date[1], cex.axis = c(0.7, 0.85))
  axis(side = 2, mgp = c(3, 0.7, 0), las = 1, cex.axis = 0.85)
  if (length(x$peaks) > 0) {
    axis(side = 3, at = time_peaks, labels = x$peaks,
         tick = FALSE, mgp = c(3, mgp2_peaks, 0), gap.axis = 0,
         col.axis = col_peaks, cex.axis = 0.7)
    axis(side = 3, at = 0, labels = "peak index:",
         tick = FALSE, mgp = c(3, mgp2_peaks, 0), hadj = 1,
         col.axis = col_peaks, cex.axis = 0.7)
  }
  if (length(x$troughs) > 0) {
    axis(side = 3, at = time_troughs, labels = x$troughs,
         tick = FALSE, mgp = c(3, mgp2_troughs, 0), gap.axis = 0,
         col.axis = col_troughs, cex.axis = 0.7)
    axis(side = 3, at = 0, labels = "trough index:",
         tick = FALSE, mgp = c(3, mgp2_troughs, 0), hadj = 1,
         col.axis = col_troughs, cex.axis = 0.7)
  }
  title(xlab = xlab, line = 2)
  title(ylab = ylab, line = 3)
  title(main = sprintf("spar = %.3f", x$spar),
        line = 2.5, cex.main = 0.9)

  invisible(NULL)
}

#' @rdname plot.smooth_cases
#' @export
#' @import graphics
dline <- function(v, ...) {
  date <- get("date.smooth_cases", envir = .epigrowthfit)
  check(v,
    what = c("numeric", "Date", "character"),
    "`v` must have class \"numeric\", \"Date\", or \"character\"."
  )
  if (is.numeric(v)) {
    check(v,
      opt = seq_along(date),
      "Numeric `v` must be an integer in `seq_along(date)`."
    )
    v <- as.integer(v)
  } else if (inherits(v, "Date")) {
    check(v,
      opt = seq(min(date), max(date), by = 1),
      "Date `v` must not be earlier than `min(date)`\n",
      "or later than `max(date)`."
    )
    v <- sapply(v, function(x) which.min(abs(date - x)))
  } else if (is.character(v)) {
    v <- try(as.Date(v), silent = TRUE)
    check(v,
      not = "try-error",
      opt = seq(min(date), max(date), by = 1),
      "Character `v` must be coercible to Date\n",
      "between `min(date)` and `max(date).`"
    )
    v <- sapply(v, function(x) which.min(abs(date - x)))
  }

  time <- days(date[v], since = date[1])
  dots <- list(...)
  if (is.null(dots$col)) {
    col.axis <- rep("black", length.out = length(v))
  } else {
    col.axis <- rep(dots$col, length.out = length(v))
  }

  op <- par(get("par.smooth_cases", envir = .epigrowthfit))
  on.exit(par(op))

  abline(v = time, ...)
  for (i in seq_along(v)) {
    axis(side = 3, at = time[i], labels = v[i],
         tick = FALSE, mgp = c(3, 1.2, 0), gap.axis = 0,
         col.axis = col.axis[i], cex.axis = 0.7)
  }

  invisible(NULL)
}
