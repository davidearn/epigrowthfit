#' @import graphics
#' @importFrom stats median setNames
plot.egf <- function(x, ..., join = NULL,
                     what = c("interval", "cumulative"),
                     xty = c("date", "numeric"),
                     log = TRUE,
                     legend = FALSE,
                     tol = 0,
                     xlim = NULL,
                     ylim = NULL,
                     xlab = NULL,
                     ylab = NULL,
                     main = NULL,
                     style = get_style_default()) {
  dots <- list(...)
  if (!is.null(join)) {
    stop_if_not(
      is.character(join),
      length(join) == 1L,
      join %in% names(object$frame)[-(1:2)],
      join %in% names(dots),
      m = paste0(
        "`join` must be NULL or an element of\n",
        "`names(object$frame)[-(1:2)]` and `names(list(...))`."
      )
    )
  }
  if (length(dots) > 0L) {
    stop_if_not(
      vapply(dots, is.atomic, logical(1L)),
      lengths(dots) == 1L,
      names(dots) %in% names(object$frame)[-(1:2)],
      !duplicated(names(dots)),
      mapply("%in%", dots, lapply(object$frame[names(dots)], levels)),
      m = paste0(
        "`list(...)` must specify valid levels\n",
        "of factors in `object$frame` (one level per factor)."
      )
    )
  }
  if (length(object$frame) == 2L) {
    i <- 1L
  } else {
    d <- object$frame[-(1:2)]
    il <- !duplicated(object$index)
    if (length(dots) > 0L) {
      il <- il & Reduce("&", lapply(names(dots), function(s) d[[s]] == dots[[s]]))
    }
    nts <- sum(il)
    if (nts == 0L) {
      stop("`list(...)` must specify a nonempty interaction\n",
           "of the factors in `object$frame`.")
    } else if (nts > 1L) {
      if (!all(duplicated(d[il, ])[-1L])) {
        stop("`list(...)` must specify a unique interaction\n",
             "of the factors in `object$frame`.")
      }
    }
    i <- which(il)
  }




  if (length(dots) > 0L) {
    f <- function(s) d[[s]] %in% dots[[s]]
    d <- d[Reduce("&", lapply(names(dots), f)), ]
  }




  ## Optional graphical parameters
  dots <- list(...)

  inc <- match.arg(inc)
  xty <- match.arg(xty)
  for (s in c("log", "add", "legend")) {
    a <- get(s, inherits = FALSE)
    stop_if_not(
      is.logical(a),
      length(a) == 1L,
      !is.na(a),
      m = sprintf("`%s` must be TRUE or FALSE.", s)
    )
  }
  if (inc == "interval") {
    stop_if_not(
      is.numeric(tol),
      length(tol) == 1L,
      !is.na(tol),
      tol >= 0,
      m = "`tol` must be a non-negative number."
    )
  }
  stop_if_not(
    inherits(style, "list"),
    !is.null(names(style)),
    m = "`style` must be a named list."
  )


  if (inc == "cumulative") {
    check(x$data$cases[-1],
      no = anyNA,
      "Cannot calculate cumulative incidence due to missing values\n",
      "in `cases[-1]`."
    )
  }






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

    if (legend) {
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
      graphics::legend(x = lx, y = ly,
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
