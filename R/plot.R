#' Plot models of epidemic growth
#'
#' @description
#' Methods for plotting objects of class "egf_init" or "egf".
#'
#' @param x An "egf_init" or "egf" object.
#' @param inc One of `"cumulative"` and `"interval"`,
#'   indicating a type of incidence to plot.
#' @param xty One of `"Date"` and `"numeric"`,
#'   indicating how time is displayed on the bottom axis.
#' @param log A logical scalar. If `TRUE`, then incidence
#'   is displayed on a logarithmic scale on the left axis.
#' @param add A logical scalar. If `TRUE`, then plot elements are
#'   added to the current graphics device, which is assumed to
#'   have been started by a previous invocation of the plot method.
#' @param annotate A logical scalar. If `TRUE`, then a legend and
#'   a list of parameter estimates are displayed in the right margin.
#' @param tol A non-negative number, used only if `inc = "interval"`.
#'   `cases[i]` is highlighted according to `style$points_short` if
#'   `diff(time)[i] < (1-tol)*m` and according to `style$points_long`
#'   if `diff(time)[i] > (1+tol)*m`, where `m = median(diff(time))`.
#'   In both cases, the value of `diff(time)[i]` is printed above
#'   the point. Assign 0 to ensure that all deviations from `m`
#'   are highlighted. Assign `Inf` to disable highlighting.
#' @param style A list of lists defining the appearance of various
#'   plot elements (see Details 2).
#' @param ... Optional arguments specifying additional
#'   graphical parameters. Currently, only `xlim`, `ylim`,
#'   `xlab`, `ylab`, and `main` are implemented.
#'   See [graphics::plot()] and [graphics::par()] for a catalogue
#'   of graphical parameters.
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
#' If `xty = "date"`, then the bottom axis displays the dates specified
#' by `date`. If `xty = "numeric"`, then the bottom axis displays the
#' number of days since `date[1]`. Regardless of `xty`,
#' numeric coordinates are used, hence the left and right boundaries of
#' the plot region are defined by `range(time)`.
#'
#' The left axis measures interval or cumulative incidence (depending
#' on `inc`). When incidence is displayed on a logarithmic scale, zeros
#' are plotted as a positive number less than 1. They are therefore
#' distinguished from nonzero counts, which are always at least 1.
#'
#' Observed data, specified by `date` and either `cases` or
#' `cumsum(cases)` (depending on `inc`), are plotted as points.
#' `cases[i]` gives the number of cases observed between `date[i]`
#' and `date[i+1]`, while `cumsum(cases)[i]` gives the number
#' observed between `date[1]` and `date[i+1]`. Both are plotted
#' at `time[i+1]`.
#'
#' The fitting window, specified by indices `first` and `last` of
#' `cases`, is displayed as a shaded rectangle behind the other plot
#' elements. The left and right boundaries occur at `time[first]`
#' and `time[last+1]`. (`cases[first]` is a count from `time[first]`
#' to `time[first+1]`, and `cases[last]` is a count from `time[last]`
#' to `time[last+1]`. Hence the fitting window starts at `time[first]`
#' and ends at `time[last+1]`.)
#'
#' The incidence curve predicted by parameter estimates `theta_init`
#' or `theta_fit` (depending on the plot method) is displayed as a line
#' supported on grid points `wgrid = seq(time[first], time[last+1], by)`,
#' where `by = 1` for cumulative incidence and `by = median(diff(time))`
#' for interval incidence (ensuring that the interval incidence curve
#' has the correct scale; see below). The predicted curve is obtained
#' with `predict(x, time = wgrid)`.
#'
#' Careful interpretation of interval incidence is required if the
#' plotted time series is not equally spaced, because `cases` roughly
#' scales with `diff(time)`. That is, certain observations may vary
#' from the predicted curve not due to chance, but because they
#' represent a count over fewer or more days than the typical
#' observation interval, namely `median(diff(time))`. Observations for
#' which `diff(time)` differs from `median(diff(time))` are highlighted
#' according to argument `tol` and labeled with the value of `diff(time)`.
#'
#' If `annotate = TRUE` and `add = FALSE`, then a legend and the
#' parameter estimates `theta_init` or `theta_fit` (depending on
#' the plot method) are displayed in the right margin.
#'
#' The plot method for class "egf" displays, in addition,
#' the doubling time associated with the "egf" object and
#' the associated 95% confidence interval.
#'
#' ## 2. Customization
#'
#' The appearance of most of the plot elements can be controlled
#' using argument `style`. `style` must be a list containing some
#' subset of the elements below.
#'
#' \describe{
#'   \item{`date`}{
#'     A named list of arguments to [daxis()]
#'     (a subset of `tcl`, `mgp2`, `col.axis`, and `cex.axis`),
#'     affecting the appearance of the bottom axis if `xty = "date"`.
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
#'     are implemented. `text_dbl` can further specify coordinates
#'     `x` and `y`.
#'   }
#' }
#'
#' List elements not specified by `style` are taken from
#' `get_style_default(x)`.
#'
#' @seealso [egf_init()], [egf()]
#' @name plot.egf
NULL

#' @rdname plot.egf
#' @export
#' @importFrom graphics plot
plot.egf_init <- function(x, inc = "interval", xty = "date", log = TRUE,
                          add = FALSE, annotate = FALSE, tol = 0,
                          style = get_style_default(x), ...) {
  ## Disguise `x` as an "egf" object to reuse `plot.egf()` machinery
  x <- list(
    init = x,
    theta_fit = x$theta_init,
    eval_cum_inc = x$eval_cum_inc,
    init_flag = TRUE
  )
  x <- structure(x, class = c("egf", "list"))
  plot(x, inc = inc, xty = xty, log = log,
       add = add, annotate = annotate, tol = tol,
       style = style, ...)
}

#' @rdname plot.egf
#' @export
#' @import graphics
#' @importFrom stats median
plot.egf <- function(x, inc = "interval", xty = "date", log = TRUE,
                     add = FALSE, annotate = FALSE, tol = 0,
                     style = get_style_default(x), ...) {
  check(inc,
    what = "character",
    len = 1,
    opt = c("interval", "cumulative"),
    "`inc` must be one of \"interval\", \"cumulative\"."
  )
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
      opt = c(TRUE, FALSE),
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

  ## Optional graphical parameters
  dots <- list(...)

  ## Observed data
  data <- data.frame(
    time = x$init$time[-1],
    cum_inc = cumsum(x$init$cases),
    int_inc = x$init$cases,
    dt = diff(x$init$time)
  )
  m <- median(data$dt)

  ## Predicted curve
  wleft <- x$init$time[x$init$first]
  wright <- x$init$time[x$init$last+1]
  wgrid <- seq(wleft, wright, by = if (inc == "interval") m else 1)
  wpred <- predict(x, wgrid)[c("time", "cum_inc", "int_inc")]
  wpred$int_inc <- c(NA, wpred$int_inc)

  ## A way to mostly avoid conditional `if (inc = ...) ... else ...`
  varname <- substr(inc, start = 1, stop = 3)
  varname <- paste0(varname, "_inc")
  formula <- as.formula(paste(varname, "~ time"))

  ## A way to include zeros on a logarithmic scale
  ymax <- max(data[[varname]])
  zero <- if (log) ymax^-0.04 else 0
  data[[varname]][data[[varname]] == 0] <- zero

  ## Titles
  if ("xlab" %in% names(dots)) {
    xlab <- dots$xlab
  } else if (tolower(xty) == "date") {
    xlab <- "date"
  } else if (xty == "numeric") {
    xlab <- paste("days since", as.character(x$init$date[1]))
  }
  if ("ylab" %in% names(dots)) {
    ylab <- dots$ylab
  } else {
    ylab <- paste(inc, "incidence")
  }
  cstr <- x$init$curve
  substr(cstr, 1, 1) <- toupper(substr(cstr, 1, 1))
  if ("main" %in% names(dots)) {
    main <- dots$main
  } else {
    main <- paste(cstr, "model of", inc, "incidence\n(fitted)")
  }

  ## Axis limits (x)
  if ("xlim" %in% names(dots)) {
    xlim <- dots$xlim
  } else {
    xmin <- 0
    xmax <- max(x$init$time) * 1.04
    xlim <- c(xmin, xmax)
  }

  ## Axis limits (y)
  if ("ylim" %in% names(dots)) {
    ylim <- dots$ylim
  } else {
    ymin <- zero
    ymax <- if (log) ymax^1.04 else ymax * 1.04
    ylim <- c(ymin, ymax)
  }

  ## Style
  s <- get_style_default(x)
  for (pe in names(s)) {
    if (!pe %in% names(style)) {
      next
    }
    gp <- intersect(names(s[[pe]]), names(style[[pe]]))
    s[[pe]][gp] <- style[[pe]][gp]
  }
  style <- s

  ## Points have sub-styles (for interval incidence,
  ## we want appearance to depend on observation interval)
  dt_min <- (1 - tol) * m
  dt_max <- (1 + tol) * m
  dt_enum <- 1 + (inc == "interval") *
    (1 * (data$dt < dt_min) + 2 * (data$dt > dt_max))
  ptys <- c("main", "short", "long")
  data$pty <- ptys[dt_enum]


  ### PLOT #############################################################

  if (add) {
    sp <- get("par", envir = .egf_env)
    sp$yaxp <- NULL
    op <- par(sp)
    on.exit(par(op))
  } else {
    op <- par(
      mar = c(4, 5, 2.7, 0.5 + 6 * annotate) + 0.1,
      las = 1,
      mgp = c(3, 0.7, 0)
    )
    on.exit({
      assign("par", par(no.readonly = TRUE), envir = .egf_env)
      par(op)
    })
    plot.new()
    plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i",
                log = if (log) "y" else "")
  }

  ## Fitting window
  l <- list(x = c(wleft, wright, wright, wleft), y = ylim[c(1, 1, 2, 2)])
  do.call(polygon, c(l, style$window))
  windex <- (data$time >= wleft - 4 & data$time <= wright + 4)

  if (!add) {
    ## Box
    box(bty = "l")

    ## Axis (x)
    if (xty == "date") {
      l <- list(
        left = par("usr")[1],
        right = par("usr")[2],
        refdate = x$init$date[1]
      )
      do.call(daxis, c(l, style$date))
    } else if (xty == "numeric") {
      axis(side = 1, cex.axis = 0.85)
    }

    ## Axis (y)
    yax_at <- axTicks(side = 2)
    if (max(yax_at) < 1e05) {
      yax_labels <- TRUE
      digits <- 0
    } else {
      mp <- matrix(unlist(strsplit(sprintf("%.6e", yax_at), "e")),
                   ncol = 2, byrow = TRUE)
      digits <- max(nchar(sub("0+$", "", mp[, 1]))) - 2
      man <- sprintf(paste0("%.", digits, "e"), as.numeric(mp[, 1]))
      pow <- as.character(as.numeric(mp[, 2]))
      if (all(as.numeric(man) %in% c(0, 1))) {
        yax_labels <- parse(text = paste0("10^", pow))
      } else {
        yax_labels <- parse(text = paste0(man, " %*% 10^", pow))
      }
      if (0 %in% yax_at) {
        yax_labels[yax_at == 0] <- expression(0)
      }
    }
    axis(side = 2, at = yax_at, labels = yax_labels,
         cex.axis = if (digits > 1) 0.65 else 0.85)
  }

  ## Observed data
  for (pty in ptys) {
    l <- list(
      formula = formula,
      data = data,
      subset = (data$pty == pty) & (!add | windex),
      xpd = !any(c("xlim", "ylim") %in% names(dots))
    )
    do.call(points, c(l, style[[paste0("points_", pty)]]))
  }

  ## Annotation above exceptional points
  if (any(dt_enum != 1)) {
    l <- list(
      formula = int_inc ~ time,
      data = data,
      labels = data$dt,
      subset = (dt_enum != 1) & (!add | windex)
    )
    do.call(text, c(l, style$text_hl))
  }

  ## Predicted curve
  l <- list(formula = formula, data = wpred)
  do.call(lines, c(l, style$lines))

  ## Doubling time
  if (is.null(x$init_flag) && inc == "interval") {
    estimate <- compute_doubling_time(x)
    sink(nullfile())
    ci <- confint(x, parm = "r", level = 0.95, method = "linear")
    sink(NULL)
    ci <- unname(rev(compute_doubling_time(ci)))
    dblstr <- sprintf("doubling time:\n%.1f (%.1f, %.1f) days",
                      estimate, ci[1], ci[2])
    wrange <- range(wpred$int_inc, na.rm = TRUE)
    if (is.na(style$text_dbl$y)) {
      if (log) {
        style$text_dbl$y <- wrange[1] * 10^(0.25 * diff(log10(wrange)))
      } else {
        style$text_dbl$y <- wrange[1] + 0.25 * diff(wrange)
      }
    }
    if (is.na(style$text_dbl$x)) {
      style$text_dbl$x <- min(wpred$time[wpred$int_inc > style$text_dbl$y], na.rm = TRUE)
    }
    l <- list(labels = dblstr, xpd = NA)
    do.call(text, c(l, style$text_dbl))
  }

  if (!add) {
    ## Titles
    title(xlab = xlab, line = 3)
    title(ylab = ylab, line = 4)
    title(main = main, line = 1, cex.main = 0.9)

    if (annotate) {
      ## Parameter estimates
      pstr1 <- paste0(names(x$theta_fit), " = ")
      pstr2 <- round(x$theta_fit, digits = 4)
      if ("K" %in% names(x$theta_fit)) {
        pstr2[["K"]] <- round(pstr2[["K"]])
      }
      px <- par("usr")[2] + 0.02 * diff(par("usr")[1:2])
      px <- px + max(strwidth(pstr1, cex = 0.7))
      py <- par("usr")[3] + 0.02 * diff(par("usr")[3:4])
      if (log) {
        py <- 10^py
      }
      text(px, py, paste(pstr1, collapse = "\n"),
           adj = c(1, 0), xpd = NA, cex = 0.7)
      text(px, py, paste(pstr2, collapse = "\n"),
           adj = c(0, 0), xpd = NA, cex = 0.7)

      ## Legend (beware: some ugly hacks here)
      lx <- par("usr")[2] + 0.02 * diff(par("usr")[1:2])
      ly <- par("usr")[4] - 0.02 * diff(par("usr")[3:4])
      if (log) {
        ly <- 10^ly
      }
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
             pch = c(style$points_main$pch, style$points_short$pch, style$points_long$pch, NA)[index],
             pt.bg = c(style$points_main$bg, style$points_short$bg, style$points_long$bg, NA)[index],
             lty = c(NA, NA, NA, style$lines$lty)[index],
             lwd = c(NA, NA, NA, style$lines$lwd)[index],
             col = c(style$points_main$col, style$points_short$col, style$points_long$col, style$lines$col)[index])
    }
  }

  invisible(NULL)
}

#' @rdname plot.egf
#' @keywords internal
#' @export
get_style_default <- function(x) {
  if (inherits(x, c("egf_init", "egf"))) {
    list(
      date = list(tcl = -0.2, mgp2 = c(0.05, 1), col.axis = c("black", "black"), cex.axis = c(0.7, 0.85)),
      points_main = list(pch = 21, col = "#BBBBBB", bg = "#DDDDDD", cex = 1),
      points_short = list(pch = 1, col = "#882255", bg = NA, cex = 1),
      points_long = list(pch = 16, col = "#882255", bg = NA, cex = 1),
      lines = list(lty = 1, lwd = 3, col = "#44AA99"),
      window = list(col = "#DDCC7740", border = NA),
      confband = list(col = "#44AA9940", border = NA),
      text_hl = list(pos = 3, offset = 0.3, col = "#BBBBBB", cex = 0.7, font = 2),
      text_dbl = list(x = NA, y = NA, pos = 4, offset = 0.8, col = "black", cex = 0.7, font = 1)
    )
  } else {
    stop("Class not implemented.")
  }
}

#' Plot simulations
#'
#' @description
#' A method for plotting simulated incidence curves specified
#' by objects of class "egf_sim".
#'
#' @param x An "egf_init" or "egf" object.
#' @param inc One of `"cumulative"` and `"interval"`,
#'   indicating a type of incidence to plot.
#' @param col_sim,col_pred Character or numeric scalars specifying
#'   colours for simulated and predicted incidence curves, respectively.
#' @param ... Optional arguments specifying additional
#'   graphical parameters. Currently, only `xlim`, `ylim`,
#'   `xlab`, `ylab`, and `main` are implemented.
#'   See [`plot()`][graphics::plot()] and [`par()`][graphics::par()]
#'   for a catalogue of graphical parameters.
#'
#' @return
#' `NULL` (invisibly).
#'
#' @details
#' ## Plot elements
#'
#' The bottom axis measures the number of days since
#' `x$object$init$date[1]`. The left axis measures interval
#' or cumulative incidence (depending on `inc`). The predicted
#' incidence curve is obtained as `pred$cum_inc` or `pred$int_inc`,
#' where `pred = predict(x$object, time = x$time)`. The simulated
#' incidence curves are obtained as the columns of `x$cum_inc` or
#' `x$int_inc`. These are plotted together behind the predicted
#' incidence curve. To help visualize the distribution of simulated
#' incidence at a given time, assign `col_sim` a sufficiently
#' transparent colour.
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

  ## Simulated data
  x$int_inc <- rbind(NA, x$int_inc)

  ## Predicted curve
  pred <- predict(x$object, time = x$time)[c("time", "cum_inc", "int_inc")]
  pred$int_inc <- c(NA, pred$int_inc)

  ## A way to avoid conditional `if (inc = ...) ... else ...`
  varname <- substr(inc, start = 1, stop = 3) # first three characters
  varname <- paste0(varname, "_inc")

  ## Titles
  if ("xlab" %in% names(dots)) {
    xlab <- dots$xlab
  } else {
    xlab <- paste("days since", as.character(x$object$init$date[1]))
  }
  if ("ylab" %in% names(dots)) {
    ylab <- dots$ylab
  } else {
    ylab <- paste(inc, "incidence")
  }
  cstr <- x$object$init$curve
  substr(cstr, 1, 1) <- toupper(substr(cstr, 1, 1)) # capitalize first letter
  if ("main" %in% names(dots)) {
    main <- dots$main
  } else {
    main <- paste0(cstr, " model of ", inc, " incidence\n",
                   "(", ncol(x[[varname]]), " simulations)")
  }

  ## Axis limits (x)
  if ("xlim" %in% names(dots)) {
    xlim <- dots$xlim
  } else {
    xmin <- 0
    xmax <- max(x$time) * 1.04
    xlim <- c(xmin, xmax)
  }

  ## Axis limits (y)
  if ("ylim" %in% names(dots)) {
    ylim <- dots$ylim
  } else {
    ymin <- 0
    ymax <- max(x[[varname]], na.rm = TRUE) * 1.04
    ylim <- c(ymin, ymax)
  }


  ### PLOT #############################################################

  op <- par(
    mar = c(4, 5, 2.7, 0.5) + 0.1,
    las = 1,
    mgp = c(3, 0.7, 0)
  )
  plot.new()
  plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i")

  ## Simulated data
  for (j in 1:ncol(x[[varname]])) {
    lines(x$time, x[[varname]][, j], lwd = 2, col = col_sim)
  }

  ## Predicted curves
  lines(pred$time, pred[[varname]], lwd = 3, col = col_pred)

  ## Box
  box(bty = "l")

  ## Axis (x)
  axis(side = 1, cex.axis = 0.85)

  ## Axis (y)
  yax_at <- axTicks(side = 2)
  if (max(yax_at) < 1e05) {
    yax_labels <- TRUE
    digits <- 0
  } else {
    mp <- matrix(unlist(strsplit(sprintf("%.6e", yax_at), "e")),
                 ncol = 2, byrow = TRUE)
    digits <- max(nchar(sub("0+$", "", mp[, 1]))) - 2
    man <- sprintf(paste0("%.", digits, "e"), as.numeric(mp[, 1]))
    pow <- as.character(as.numeric(mp[, 2]))
    if (all(as.numeric(man) %in% c(0, 1))) {
      yax_labels <- parse(text = paste0("10^", pow))
    } else {
      yax_labels <- parse(text = paste0(man, " %*% 10^", pow))
    }
    if (0 %in% yax_at) {
      yax_labels[yax_at == 0] <- expression(0)
    }
  }
  axis(side = 2, at = yax_at, labels = yax_labels,
       cex.axis = if (digits > 1) 0.65 else 0.85)

  ## Titles
  title(xlab = xlab, line = 3)
  title(ylab = ylab, line = 4)
  title(main = main, line = 1, cex.main = 0.9)

  par(op)
  invisible(NULL)
}

#' Plot smoothed incidence data
#'
#' @description
#' A method for plotting objects of class "smooth_cases".
#'
#' @param x A "smooth_cases" object.
#' @param ... Unused optional arguments.
#'
#' @return
#' `NULL` (invisibly).
#'
#' @details
#' ## Plot elements
#' There is one panel for each element of `x$spar`. The bottom axis
#' displays the dates specified by `x$date`. The left axis measures
#' interval incidence, log(1+x)-transformed if `x$log = TRUE`.
#' Plotted as points are the interval incidence data specified by
#' `x$date` and `x$cases`. Plotted as a solid line is a cubic spline
#' fit to the (possibly transformed) data using [stats::smooth.spline()].
#' The spline in panel `i` is specified by `x$ss[[i]]` and uses
#' smoothing parameter `x$spar[i]`. Dotted lines are drawn at peaks
#' (red) and troughs (blue) in the spline. These lines are labeled in
#' the top margin with the index of `cases` corresponding to the time
#' of the peak or trough in the spline.
#'
#' @export
#' @import graphics
plot.smooth_cases <- function(x, ...) {
  object <- x
  op <- par(
    mfrow = if (length(object$spar) < 4) c(1, 1) else c(2, 2),
    mar = c(3, 4, 3, 1) + 0.5
  )
  on.exit(par(op))

  data <- data.frame(
    x = object$time[-1],
    y = if (object$log) log10(1 + object$cases) else object$cases
  )
  xlim <- range(object$time)
  ylim <- c(0, max(data$y) * 1.04)
  xlab <- "date"
  ylab <- if (object$log) expression(log[10] * "(1+cases)") else "cases"
  col_points <- "#CCCCCC"
  col_peaks <- "#BB5566"
  col_troughs <- "#004488"

  for (i in seq_along(object$spar)) {
    peaks <- object$peaks[[i]]
    troughs <- object$troughs[[i]]
    x_peaks <- data$x[peaks]
    x_troughs <- data$x[troughs]
    mgp2_peaks <- if (length(troughs) > 0) 0.6 else 0
    mgp2_troughs <- 0

    plot.new()
    plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i")
    points(y ~ x, data = data, col = col_points)
    lines(y ~ x, data = predict(object$ss[[i]]), lwd = 2)
    abline(v = x_peaks, lty = 3, col = col_peaks)
    abline(v = x_troughs, lty = 3, col = col_troughs)
    box(bty = "l")
    l <- list(
      left = par("usr")[1],
      right = par("usr")[2],
      refdate = x$init$date[1]
    )
    daxis(left = par("usr")[1], right = par("usr")[2],
          refdate = object$date[1], cex.axis = c(0.7, 0.85))
    axis(side = 2, mgp = c(3, 0.7, 0), las = 1, cex.axis = 0.85)
    if (length(x_peaks) > 0) {
      axis(side = 3, at = x_peaks, labels = peaks,
           tick = FALSE, mgp = c(3, mgp2_peaks, 0), gap.axis = 0,
           col.axis = col_peaks, cex.axis = 0.7)
      axis(side = 3, at = 0, labels = "peak index:",
           tick = FALSE, mgp = c(3, mgp2_peaks, 0), hadj = 1,
           col.axis = col_peaks, cex.axis = 0.7)
    }
    if (length(x_troughs) > 0) {
      axis(side = 3, at = x_troughs, labels = troughs,
           tick = FALSE, mgp = c(3, mgp2_troughs, 0), gap.axis = 0,
           col.axis = col_troughs, cex.axis = 0.7)
      axis(side = 3, at = 0, labels = "trough index:",
           tick = FALSE, mgp = c(3, mgp2_troughs, 0), hadj = 1,
           col.axis = col_troughs, cex.axis = 0.7)
    }
    title(xlab = xlab, line = 2)
    title(ylab = ylab, line = 2)
    title(main = sprintf("spar = %.3f", object$spar[i]),
          line = 2, cex.main = 0.9)
  }

  invisible(NULL)
}

