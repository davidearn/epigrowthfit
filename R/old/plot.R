#' Plot models of epidemic growth
#'
#' @description
#' Methods for plotting `"egf"` objects.
#'
#' @param x
#'   An `"egf"` object returned by [egf()].
#' @param type
#'   A character string indicating a type of plot. The options are:
#'   interval incidence (`"interval"`),
#'   cumulative incidence (`"cumulative"`),
#'   per capita growth rate as curve (`"rt1"`), and
#'   per capita growth rate as heat map (`"rt2"`).
#' @param subset
#'   An expression to be evaluated in the combined model frame.
#'   Must evaluate to a logical vector or list of logical vectors
#'   indexing rows of the model frame. Only indexed fitting windows
#'   are plotted. The default (`NULL`) is to plot all fitting windows.
#' @param order
#'   An expression to be evaluated in the combined model frame,
#'   typically a call to [order()], determining the order in which
#'   time series are plotted. Must evaluate to a permutation of
#'   `seq_len(nrow)`. The default (`NULL`) is equivalent to
#'   `seq_len(nrow)`.
#' @param cache
#'   An `"egf_plot_cache"` object returned by a previous evaluation of
#'   `plot.egf(x)`. Predicted values and confidence intervals stored
#'   in `cache` will be reused (rather than recomputed) if relevant to
#'   the current function call, i.e., to `type` and `subset`.
#' @param time_as
#'   A character string indicating how time is displayed on the
#'   bottom axis. The options are: as a calendar (`"Date"`) and
#'   as a number of days since the earliest time point (`"numeric"`).
#' @param log
#'   A logical scalar. If `TRUE`, then the dependent variable is
#'   plotted on a logarithmic scale. Unused by `type = "rt[12]"`.
#' @param level
#'   A number in the interval (0,1). The desired confidence level
#'   when `show_tdoubling = TRUE` or `show_bands = TRUE`.
#' @param show_tdoubling
#'   A logical scalar. If `TRUE`, then confidence intervals on initial
#'   doubling times
#' @param show_bands
#'   A logical scalar. If `TRUE`, then confidence bands on predicted
#'   curves are drawn. Unused by `type = "rt2"`.
#' @param show_legend
#'   A logical scalar. If `TRUE`, then a legend is displayed in the
#'   right margin when `type = "interval"`.
#' @param per_plot
#'   A positive integer. The number of panels displayed in one plot
#'   when `type = "rt2"`.
#' @param control
#'   A named list controlling the appearance of various plot elements,
#'   or otherwise `NULL` (see Details).
#' @param ...
#'   Optional arguments specifying additional graphical parameters
#'   to be recycled for all plots. Currently, only `xlim`, `ylim`,
#'   `xlab`, `ylab`, and `main` are used (see Details).
#'
#' @return
#' An `"egf_plot_cache"` object. If `cache` was used, then this object
#' is the result of augmenting `cache` with any new computations.
#'
#' @details
#' `plot.egf()` will _not_ detect mismatch between `x` and `cache`.
#' In other words, constructions like `plot(x2, cache = plot(x1))`
#' should _not_ be expected to produce correct results.
#'
#' `control` must be `NULL` or a named list defining some subset of
#' the available control parameters (see section "Control parameters").
#' A list of modifiable options for each control parameter and their
#' default values can be obtained with
#' `get_control_default("plot.egf", type = type)`. Unspecified options
#' take their values in this list. Unsupported options are silently
#' discarded. To suppress a plot element, set the corresponding control
#' parameter to `NULL`, as in `control = list(box = NULL)`.
#'
#' If `xty = "Date"`, then `xlim` can be a Date vector or character
#' vector coercible to Date via `as.Date(xlim)`.
#'
#' `main` is evaluated as an expression in the combined model frame,
#' so it can be used to define a template for all plot titles. For
#' example, if plotted time series correspond to levels of a factor
#' named `country` found in the combined model frame, then one might
#' set `main = sprintf("disease incidence in %s", country)` or simply
#' `main = country`.
#'
#' Characters appearing before and after the first instance of `"\n"`
#' in `main` are displayed according to control parameters `main` and
#' `sub`, respectively (see section "Control parameters").
#'
#' @inheritSection fitted.egf Warning
#'
#' # Control parameters
#'
#' For `type != "rt2"`:
#' \describe{
#' \item{`box`}{
#'   A named list of arguments to [graphics::box()], affecting
#'   the appearance of the box drawn around the plot region.
#' }
#' \item{`xax_*`, `yax`}{
#'   Named lists of arguments to [graphics::axis()], affecting
#'   the appearance of the bottom and left axes.
#' }
#' \item{`main`, `sub`, `xlab`, `ylab`}{
#'   Named lists of arguments to [graphics::title()], affecting
#'   the appearance of axis titles.
#' }
#' \item{`rect`}{
#'   A named list of arguments to [graphics::rect()], affecting
#'   the appearance of fitting windows.
#' }
#' \item{`polygon`}{
#'   A named list of arguments to [graphics::polygon()], affecting
#'   the appearance of confidence bands on predicted curves.
#' }
#' \item{`lines`}{
#'   A named list of arguments to [graphics::lines()], affecting
#'   the appearance of predicted curves.
#' }
#' \item{`points`}{
#'   A named list of arguments to [graphics::points()], affecting
#'   the appearance of observed data.
#' }
#' \item{`text_tdoubling_*`}{
#'   Named lists of arguments to [graphics::text()], affecting
#'   the appearance of initial doubling times printed when
#'   `x$curve %in% c("exponential", "logistic", "richards")`.
#'   The caption, estimate, and confidence interval are controlled
#'   separately.
#' }
#' }
#' For `type = "interval"`:
#' \describe{
#' \item{`tol`}{
#'   A non-negative number or `Inf`, defining "exceptional" points.
#'   Within a time series, `x[-1]` is highlighted
#'   according to `points_short` if `diff(time) < (1-tol)*m` and
#'   according to `points_long`  if `diff(time) > (1+tol)*m`,
#'   where `m = median(diff(time))`. In both cases, the value
#'   of `diff(time)` is printed above the point according to
#'   `text_short_long`.
#'   Assign 0 to highlight all deviations from `m`.
#'   Assign `Inf` to disable highlighting.
#' }
#' \item{`points_short`,`points_long`}{
#'   Alternatives to `points` used for exceptional points. See `tol`.
#' }
#' \item{`text_short_long`}{
#'   A named list of arguments to [graphics::text()], affecting
#'   the appearance of text above exceptional points. See `tol`.
#' }
#' }
#' For `type = "rt1"`:
#' \describe{
#' \item{`abline`}{
#'   A named list of arguments to [graphics::abline()],
#'   affecting the appearance of the line drawn at `y = 0`.
#' }
#' \item{`segments`}{
#'   A named list of arguments to [graphics::segments()],
#'   affecting the appearance of line segments drawn at `y = r`
#'   when `x$curve %in% c("exponential", "logistic", "richards")`.
#' }
#' }
#' For `type = "rt2"`:
#' \describe{
#' \item{`colorRamp`}{
#'   A named list of arguments to [grDevices::colorRamp()],
#'   defining the heat map's color palette.
#' }
#' \item{`ips`}{
#'   A non-negative number, defining the space between panels
#'   as a number of margin lines.
#' }
#' }
#'
#' @export
plot.egf <- function(x,
                     type = c("interval", "cumulative", "rt1", "rt2"),
                     subset = NULL,
                     order = NULL,
                     cache = NULL,
                     time_as = c("Date", "numeric"),
                     log = TRUE,
                     level = 0.95,
                     show_tdoubling = FALSE,
                     show_bands = FALSE,
                     show_legend = FALSE,
                     per_plot = 6L,
                     control = NULL,
                     ...) {
  type <- match.arg(type)
  time_as <- match.arg(time_as)
  control <- get_clean_control(control, "plot.egf", type = type)
  if (!grepl("^rt[12]$", type)) {
    stop_if_not_true_false(log)
  }
  if (type == "rt2") {
    stop_if_not_positive_integer(per_plot)
  } else {
    if (x$curve %in% c("exponential", "logistic", "richards") &&
        (!is.null(control$text_tdoubling_estimate) || !is.null(control$text_tdoubling_ci))) {
      stop_if_not_true_false(show_tdoubling)
    } else {
      show_tdoubling <- FALSE
    }
    if (!is.null(control$bands)) {
      stop_if_not_true_false(show_bands)
    } else {
      show_bands <- FALSE
    }
    if (type == "interval") {
      stop_if_not_true_false(show_legend)
    }
  }

  frame <- do.call(cbind, unname(object$frame_par))
  frame <- frame[!duplicated(names(frame))]
  subset <- subset_bak <- subset_to_index(substitute(subset), frame, parent.frame())
  order <- order_to_index(substitute(order), frame, parent.frame())

  what <- switch(type, interval = "log_int_inc", cumulative = "log_cum_inc", "log_rt")
  window <- x$frame_ts$window
  need <- levels(window)[subset]
  if (!is.null(cache)) {
    stop_if_not(
      inherits(cache, "egf_plot_cache"),
      m = "`cache` must be NULL or inherit from class \"egf_plot_cache\"."
    )
    useful <- (cache$var == sub("^log_", "", what)) & (!show_bands | !is.na(cache$lower))
    have <- levels(droplevels(cache$window[useful]))
    need <- setdiff(need, have)
    subset <- match(need, levels(window), 0L)
  }
  cn <- c("var", "ts", "window", "time" "estimate", "lower", "upper")
  if (length(need) > 0L) {
    window <- factor(window, levels = need)
    if (type == "interval") {
      dt <- tapply(time, window, function(x) median(diff(x)), 0)
    } else {
      dt <- 1
    }
    time_split <- Map(seq.int,
      from = x$endpoints$start[subset],
      to = x$endpoints$end[subset],
      by = dt
    )
    pd <- predict(x,
      what = what,
      time = unlist(time_split),
      window = rep.int(x$endpoints$window[subset], lengths(time_split)),
      log = show_bands,
      se = show_bands
    )
    if (show_bands) {
      cache <- rbind(cache, confint(pd, level = level, log = FALSE)[cn])
    } else {
      pd$lower <- pd$upper <- NA_real_
      cache <- rbind(cache, pd[cn])
    }
    cache <- cache[do.call(order, cache[cn]), , drop = FALSE]
    cache <- cache[!duplicated(cache[cn[1:5]]), , drop = FALSE]
  }
  if (show_tdoubling && !"tdoubling" %in% levels(cache$var)) {
    ci <- confint(x, level = level, par = "tdoubling", method = "wald")
    names(ci) <- sub("^par$", "var", names(ci))
    ci$time <- NA_real_
    cache <- rbind(ci[cn], cache)
  }

  dots <- list(...)


  if (type == "rt2") {
    plot.egf_heat(x,
      cache = cache,
      time_as = time_as,
      per_plot = per_plot,
      control = control,
      ...
    )
  } else {
    plot.egf_curve(x,
      type = type,
      subset = subset_bak,
      order = order,
      cache = cache,
      time_as = time_as,
      log = log,
      level = level,
      show_tdoubling = show_tdoubling,
      show_bands = show_bands,
      show_legend = show_legend,
      control = control,
      ...
    )
  }

  class(cache) <- c("egf_plot_cache", "data.frame")
  cache
}

#' @import graphics
#' @importFrom stats reformulate median predict confint
plot.egf.main <- function(x, type, log, tol, legend, bands, level,
                          control, frame_aug_split, ...) {
  ### Validation ==========================================

  dots <- list(...)
  if (type == "cumulative") {
    omit <- vapply(frame_aug_split, function(d) anyNA(d[-1L, 2L]), FALSE)
    if (any(omit)) {
      warning("Missing values in interval incidence time series\n",
              "preventing calculation of cumulative incidence.")
    }
    if (all(omit)) {
      stop("There were no time series without missing values.")
    }
    frame_aug_split <- frame_aug_split[!omit]
  }
  if (type == "interval") {
    stop_if_not(
      is.numeric(tol),
      length(tol) == 1L,
      tol >= 0,
      m = "`tol` must be a non-negative number or `Inf`."
    )
  }

  ### Set-up ==============================================

  ## A way to avoid conditional `if (type == ...) expr1 else expr2`
  ## in some places
  varname <- switch(type,
    interval = "int_inc",
    cumulative = "cum_inc",
    rt1 = "rt"
  )
  formula <- reformulate("time", varname)

  ## Confidence intervals on initial doubling times
  show_tdoubling <- (
    type %in% c("interval", "rt1") &&
    x$curve %in% c("exponential", "logistic", "richards")
  )
  if (show_tdoubling) {
    ci <- confint(x, parm = "tdoubling", level = level, method = "wald")
  }

  ## Axis title (y)
  if (is.null(dots$ylab)) {
    if (type %in% c("interval", "cumulative")) {
      ylab <- sprintf("%s incidence", type)
    } else {
      ylab <- "growth rate, per day"
      ylab_again <- "doubling time, days"
    }
  } else {
    if (type %in% c("interval", "cumulative")) {
      ylab <- dots$ylab
    } else {
      ylab <- dots$ylab[1L]
      ylab_again <- dots$ylab[2L]
    }
  }

  ## Axis title (main)
  group_by <- attr(x$frame, "group_by")
  if (is.null(dots$main)) {
    if (type %in% c("interval", "cumulative")) {
      s <- switch(x$curve, gompertz = "Gompertz", richards = "Richards", x$curve)
      main_template <- sprintf("Fitted %s model", s)
    } else {
      main_template <- "Instantaneous exponential growth rate"
    }
    if (length(group_by) > 0L) {
      s <- paste(sprintf("%s = %%%s", group_by, group_by), collapse = ", ")
      main_template <- sprintf("%s\n%s", main_template, s)
    }
  } else {
    main_template <- dots$main
  }

  ### Loop over plots =====================================

  mar <- switch(type,
    rt1 = c(3, 4 + (is.null(dots$ylim) || dots$ylim[2L] > 0) * 4, 4, 1) + 0.1,
    c(3, 5, 4, 1 + 5.5 * legend) + 0.1
  )
  op <- par(mar = mar)
  on.exit(par(op))

  for (d in frame_aug_split) {

    ### Set-up for plot ===================================

    index <- d$.index
    date <- d[[1L]]
    cases <- d[[2L]]
    i12 <- vapply(levels(index), function(s) range(which(index == s)), integer(2L))
    i1 <- i12[1L, ]
    i2 <- i12[2L, ]
    d0 <- date[1L]
    d1 <- date[i1]
    d2 <- date[i2]
    t1 <- days(d1, since = d0)
    t2 <- days(d2, since = d0)

    ## Data frame for plot
    data <- data.frame(
      time = days(date, since = d0),
      cum_inc = cumsum(c(0L, cases[-1L])),
      int_inc = cases
    )
    data$dt <- c(NA, diff(data$time))
    data$rt <- c(diff(base::log(data$int_inc)) / data$dt[-1L], NA)
    data$rt[!is.finite(data$rt)] <- NA

    ## Predicted curves with confidence bands
    m <- median(data$dt[-1L])
    make_pred <- function(i1, t1, t2, by) {
      p <- predict(x,
        time = seq.int(from = 0, to = t2 - t1, by = by),
        varname = sprintf("log_%s", varname),
        subset = if (length(d) > 3L) d[i1, -c(1:2, length(d)), drop = FALSE],
        se = bands
      )
      if (bands) {
        p <- confint(p, level = level, log = TRUE)
      }
      p <- p[[1L]]
      p[[1L]] <- t1 + p[[1L]]
      p[-1L] <- exp(p[-1L])
      if (type == "cumulative" && i1 > 1L) {
        p[-1L] <- sum(cases[2L:i1]) + p[-1L]
      }
      p
    }
    pred <- Map(make_pred, i1 = i1, t1 = t1, t2 = t2,
                by = switch(type, interval = m, 1))

    ## Title with flags substituted
    if (is.character(main_template)) {
      ## Replace "%factor_name" with "level_name"
      for (s in group_by) {
        main <- gsub(sprintf("%%%s", s), as.character(d[1L, s]), main_template, fixed = TRUE)
      }
    } else {
      main <- main_template
    }

    ## Axis limits (x)
    if (is.null(dots$xlim)) {
      xlim <- c(0, max(data$time) * 1.04)
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
      if (type %in% c("interval", "cumulative")) {
        ymax <- max(data[[varname]], na.rm = TRUE)
        if (log) {
          zero <- ymax^-0.04
          data[[varname]][data[[varname]] == 0L] <- zero
        } else {
          zero <- 0
        }
        ylim <- c(zero, if (log) ymax^1.04 else ymax * 1.04)
      } else {
        ylim <- range(data$rt, na.rm = TRUE)
        ylim[1L] <- max(-log(2), ylim[1L])
        ylim[2L] <- min(log(2), ylim[2L])
        ylim <- ylim + c(-1, 1) * 0.04 * diff(ylim)
      }
    } else {
      ylim <- dots$ylim
    }

    ## Point styles depending on observation interval
    control$points_median <- control$points
    pty <- c("median", "short", "long")
    if (type == "interval") {
      dt_min <- (1 - tol) * m
      dt_max <- (1 + tol) * m
      dt_enum <- 1L + 1L * (data$dt < dt_min) + 2L * (data$dt > dt_max)
      data$pty <- factor(pty)[dt_enum]
    } else {
      data$pty <- factor("median")
    }

    ### Plot ==============================================

    plot.new()
    plot.window(
      xlim = xlim,
      ylim = ylim,
      xaxs = "i",
      yaxs = "i",
      log = if (log) "y" else ""
    )
    cex <- par("cex")
    csi <- par("csi")
    cxy <- par("cxy")
    pin <- par("pin")
    usr <- par("usr")

    ## Fitting windows
    if (!is.null(control$rect)) {
      l <- list(
        xleft = t1,
        ybottom = inv_log10(usr[3L], log),
        xright = t2,
        ytop = inv_log10(usr[4L], log)
      )
      do.call(rect, c(l, control$rect))
    }

    ## Observed data
    for (s in pty) {
      if (!is.null(control[[sprintf("points_%s", s)]])) {
        l <- list(
          formula = formula,
          data = data,
          subset = (data$pty == s)
        )
        do.call(points, c(l, control[[sprintf("points_%s", s)]]))
      }
    }

    ## Annotation above exceptional points
    if (type == "interval" && !is.null(control$text_short_long)) {
      l <- list(
        formula = formula,
        data = data,
        labels = data$dt,
        subset = (data$pty != pty[1L])
      )
      do.call(text, c(l, control$text_short_long))
    }

    ## Confidence band
    if (bands && !is.null(control$polygon)) {
      for (p in pred) {
        l <- list(
          x = c(p$time, rev(p$time)),
          y = c(p$lower, rev(p$upper))
        )
        do.call(polygon, c(l, control$polygon))
      }
    }

    ## Predicted curve
    if (!is.null(control$lines)) {
      for (p in pred) {
        l <- list(formula = estimate ~ time, data = p)
        do.call(lines, c(l, control$lines))
      }
    }

    ## Asymptotes
    if (type == "rt1") {
      if (!is.null(control$abline)) {
        do.call(abline, c(list(h = 0), control$abline))
      }
      if (!is.null(control$segments)) {
        for (s in levels(index)) {
          p <- pred[[s]]
          r <- log(2) / ci[match(s, levels(x$index)), length(ci) - 2L]
          l <- list(
            x0 = p$time[1L],
            y0 = r,
            x1 = p$time[nrow(p)],
            y1 = r
          )
          do.call(segments, c(l, control$segments))
        }
      }
    }

    ## Initial doubling times
    if (show_tdoubling) {
      control$text_tdoubling <- lapply(control$text_tdoubling, rep_len, 3L)
      l <- lapply(1:3, function(i) lapply(control$text_tdoubling, `[`, i))
      names(l) <- c("e", "ci", "cap")
      ## Much ado about finding the right user coordinates
      ## when log scale is in effect
      y_ci <- add_lines_to_user(0.25, inv_log10(usr[4L], log), log)
      h_ci <- strheight("", cex = l$ci$cex, font = l$ci$font)
      y_e <- add_lines_to_user(0.15, add_height_to_user(h_ci, y_ci, log), log)
      h_e <- strheight("", cex = l$e$cex, font = l$e$font)
      for (s in levels(index)) {
        p <- pred[[s]]
        elu <- unlist(ci[match(s, levels(x$index)), length(ci) - 2:0])
        text(
          x = rep.int(mean(p$time[c(1L, nrow(p))]), 2L),
          y = c(y_e, y_ci),
          labels = c(sprintf("%.1f", elu[1L]), sprintf("(%.1f, %.1f)", elu[2L], elu[3L])),
          cex = cex * c(l$e$cex, l$ci$cex),
          col = c(l$e$col, l$ci$col),
          font = c(l$e$font, l$ci$font),
          adj = c(0.5, 0),
          xpd = NA
        )
      }
      y_ci_again <- add_lines_to_user(0.5, add_height_to_user(h_e, y_e, log), log)
      y_e_again <- add_lines_to_user(0.15, add_height_to_user(h_ci, y_ci_again, log), log)
      text(
        x = rep.int(usr[2L] - 0.5 * strwidth("estimate", cex = control$text_tdoubling$cex, font = 2), 2L),
        y = c(y_e_again, y_ci_again),
        labels = c("estimate", "(95% CI)"),
        cex = cex * c(l$e$cex, l$ci$cex),
        col = c(l$e$col, l$ci$col),
        font = c(l$e$font, l$ci$font),
        adj = c(0.5, 0),
        xpd = NA
      )
      text(
        x = usr[2L],
        y = add_lines_to_user(0.25, add_height_to_user(h_e, y_e_again, log), log),
        labels = "initial doubling time, days:",
        cex = cex * l$cap$cex,
        col = l$cap$col,
        font = l$cap$font,
        adj = c(1, 0),
        xpd = NA
      )
    }

    ## Box
    if (!is.null(control$box)) {
      do.call(box, control$box)
    }

    ## Axis (x)
    if (!is.null(control$xax)) {
      l <- list(
        left = usr[1L],
        right = usr[2L],
        refdate = d0
      )
      do.call(daxis, c(l, control$xax))
    }

    ## Axis (y)
    if (!is.null(control$yax)) {
      l1 <- list(side = 2, las = 1)
      l2 <- control$yax
      l2$mgp <- c(3, control$yax$mgp2, 0)
      l2$mgp2 <- NULL
      l1$at <- axTicks(side = 2)
      if (type %in% c("interval", "cumulative") && (max(l1$at) >= 1e05)) {
        l1$labels <- get_labels(l1$at)
        l2$cex.axis <- min(l2$cex.axis, get_cex_axis(l1$labels, mex = 3.5 - control$yax$mgp2))
      }
      do.call(axis, c(l1, l2))
      if (type == "rt1" && usr[4L] > 0) {
        axis(
          side = 2,
          line = 4.5,
          at = c(0, usr[4L]),
          labels = c("", ""),
          lwd.ticks = 0
        )
        tdoubling <- c(1:5, 10, 50, 100)
        l1$at <- log(2) / tdoubling
        l1$labels <- tdoubling
        ## FIXME: why mgp[4] = 4 + x, not mgp[2] = 4.5 + x?
        l2$mgp <- c(3, 4 + control$yax$mgp2, 4.5)
        l2$lwd <- 0
        l2$lwd.ticks <- 1
        do.call(axis, c(l1, l2))
      }

      ## Axis title (y)
      if (!is.null(control$ylab)) {
        if (type %in% c("interval", "cumulative")) {
          line0 <- 4
        } else {
          line0 <- 0.5 + max(strwidth(axTicks(side = 2), units = "inches", cex = control$yax$cex.axis, font = control$yax$font.axis)) / csi + control$yax$mgp2
        }
        l <- list(ylab = ylab, line = line0)
        do.call(title, c(l, control$ylab))
        if (type == "rt1" && usr[4L] > 0) {
          ## String width relative to 0.8 times the distance from 0 to usr[4]
          rsw <- (strwidth(ylab_again, units = "inches", cex = control$ylab$cex.lab, font = control$ylab$font.lab) * diff(usr[3:4]) / pin[2L]) / (0.8 * (usr[4L] - max(0, usr[3L])))
          text(
            ## FIXME: why 7.5, not 4.5 + 0.5?
            x = usr[1L] - cxy[1L] * (7.5 + max(strwidth(tdoubling, units = "inches", cex = control$yax$cex.axis, font = control$yax$font.axis)) / csi + control$yax$mgp2),
            y = mean(c(max(0, usr[3L]), usr[4L])),
            labels = ylab_again,
            adj = c(0.5, 0),
            srt = 90,
            xpd = NA,
            col = control$ylab$col.lab,
            cex = control$ylab$cex.lab * cex / max(1, rsw),
            font = control$ylab$font.lab
          )
        }
      }
    }

    ## Axis title (main)
    if (!is.null(control$main)) {
      main_split <- strsplit(main, "\n", fixed = TRUE)[[1L]]
      control$main <- lapply(control$main, rep_len, 2L)
      if (show_tdoubling) {
        line0 <- user_to_lines(add_lines_to_user(0.5, add_height_to_user(h_e, y_e, log), log), log)
      } else {
        line0 <- 0.5
      }
      if (length(main_split) > 1L) {
        l2 <- lapply(control$main, `[`, 1L)
        l2_sub <- lapply(control$main, `[`, 2L)
        h <- strheight("", units = "inches", cex = control$main$cex.main[2L], font = control$main$font.main[2L]) / csi
        for (i in seq.int(length(main_split), 2L)) {
          l1 <- list(main = main_split[i], line = line0)
          do.call(title, c(l1, l2_sub))
          line0 <- line0 + h + 0.15
        }
        l1 <- list(main = main_split[1L], line = line0 + 0.1)
        do.call(title, c(l1, l2))
      } else {
        control$main <- lapply(control$main, `[`, 1L)
        l <- list(main = main, line = line0)
        do.call(title, c(l, control$main))
      }
    }

    ## Legend
    if (legend) {
      lx <- usr[2L] + 0.02 * diff(usr[1:2])
      ly <- inv_log10(usr[4L] - 0.02 * diff(usr[3:4]), log)
      cond <- all(data$dt[data$pty == "median"] == m, na.rm = TRUE)
      lexp <- parse(
        text = sprintf("'%s,' ~ Delta * 't %s %g day%s'",
          c("obs", "obs", "obs", "pred"),
          c((if (cond) "=" else "~"), "<", ">", "="),
          m,
          if (m > 1) "s" else ""
        )
      )
      lind <- c(pty %in% levels(data$pty), TRUE)
      get_el <- function(pe) sapply(control[sprintf("points_%s", pty)], `[[`, pe)
      graphics::legend(
        x = lx,
        y = ly,
        legend = lexp[lind],
        xpd = NA,
        bty = "n",
        cex = 0.7,
        seg.len = 1,
        pch = c(get_el("pch"), NA)[lind],
        pt.bg = c(get_el("bg"), NA)[index],
        lty = c(rep.int(NA, 3L), control$lines$lty)[lind],
        lwd = c(rep.int(NA, 3L), control$lines$lwd)[lind],
        col = c(get_el("col"), control$lines$col)[lind]
      )
    }
  }

  invisible(NULL)
}

#' @import graphics
#' @importFrom grDevices colorRamp rgb
plot.egf.heat <- function(x, per_plot, bw_panels,
                          control, frame_aug_split, ...) {
  dots <- list(...)
  stop_if_not_positive_integer(per_plot)
  stop_if_not(
    is.numeric(bw_panels),
    length(bw_panels) == 1L,
    bw_panels >= 0,
    m = "`bw_panels` must be a non-negative number."
  )

  group_by <- attr(x$frame, "group_by")
  refdate <- do.call(min, lapply(frame_aug_split, `[`, 1L, 1L))
  xlim <- c(0, days(do.call(max, lapply(frame_aug_split, function(d) d[nrow(d), 1L])), since = refdate))
  ylim <- c(0, 1)

  f <- function(d) {
    index <- d$.index
    i12 <- vapply(levels(index), function(s) range(which(index == s)), integer(2L))
    i1 <- i12[1L, ]
    i2 <- i12[2L, ]
    t1 <- days(d[i1, 1L], since = refdate)
    t2 <- days(d[i2, 1L], since = refdate)

    g <- function(i1, t1, t2) {
      p <- predict(x,
        time = seq.int(from = 0, to = t2 - t1, by = 1),
        varname = "log_rt",
        subset = if (length(d) > 3L) d[i1, -c(1:2, length(d)), drop = FALSE],
        se = FALSE
      )
      p <- p[[1L]]
      p$time <- t1 + p$time
      p$estimate <- exp(p$estimate)
      p
    }
    Map(g, i1 = i1, t1 = t1, t2 = t2)
  }
  pred <- lapply(frame_aug_split, f)
  rt_max <- do.call(max, lapply(do.call(c, pred), function(d) max(d$estimate, na.rm = TRUE)))

  ### Loop over plots =====================================

  op <- par(oma = c(2.5, 0, 3.5, 0))
  on.exit(par(op))

  pal <- do.call(colorRamp, control$colorRamp)

  j <- 0L
  while (j < length(frame_aug_split)) {

    layout(matrix(c(seq_len(per_plot), rep.int(per_plot + 1L, per_plot)), ncol = 2L), widths = c(par("din")[1L] - 1.2, 1.2))
    par(mar = c(0.5 * bw_panels, 1.5, 0.5 * bw_panels, 0.5))

    ### Loop over panels ==================================

    for (k in j + seq_len(min(per_plot, length(frame_aug_split) - j))) {
      plot.new()
      plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i")

      polygon(
        x = xlim[c(1L, 2L, 2L, 1L)],
        y = ylim[c(1L, 1L, 2L, 2L)],
        border = NA,
        col = "black"
      )

      ### Loop over fitting windows =======================

      for (p in pred[[k]]) {
        for (i in seq_len(nrow(p))) {
          polygon(
            x = (p$time[i] + c(-0.5, 0.5))[c(1L, 2L, 2L, 1L)],
            y = ylim[c(1L, 1L, 2L, 2L)],
            border = NA,
            col = rgb(pal(p$estimate[i] / rt_max), maxColorValue = 255)
          )
        }
      } # loop over fitting windows

      s <- names(frame_aug_split)[k]

      w <- strwidth(s, cex = 1, font = 2)
      h <- strheight(s, cex = 1, font = 2)
      x_pad <- diff(par("usr")[1:2]) * (0.075 * par("pin")[2L] / par("pin")[1L])
      y_pad <- 0.075 * diff(par("usr")[3:4])
      rect(
        xleft = par("usr")[1L],
        ybottom = par("usr")[4L] - h - 2 * y_pad,
        xright = par("usr")[1L] + w + 2 * x_pad,
        ytop = par("usr")[4L],
        border = NA,
        col = "#00000080"
      )
      text(
        x = par("usr")[1L] + x_pad,
        y = par("usr")[4L] - y_pad,
        labels = s,
        adj = c(0, 1),
        cex = 1,
        font = 2,
        col = "white"
      )
      if (k == j + 1L) {
        title(
          main = sprintf("Instantaneous exponential growth rate, by %s", paste(group_by, collapse = ":")),
          line = 0.5,
          adj = 0,
          cex.main = 1.3,
          xpd = NA
        )
      }
    } # loop over panels

    daxis(
      left = xlim[1L],
      right = xlim[2L],
      refdate = refdate,
      tcl = -0.2,
      mgp2 = c(0.25, 1.25),
      cex.axis = c(0.85, 1.15) + 0.1
    )
    for (l in seq_len((per_plot - (k %% per_plot)) %% per_plot)) {
      plot.new()
    }

    par(mar = c(0.5 * bw_panels, 1, 0.5 * bw_panels, 7))
    plot.new()
    dy2 <- 0.0025
    plot.window(
      xlim = c(0, 1),
      ylim = c(-dy2, 1 + dy2) * rt_max,
      xaxs = "i",
      yaxs = "i"
    )
    for (y in seq.int(0, 1, by = 2 * dy2)) {
      polygon(
        x = c(0, 1)[c(1L, 2L, 2L, 1L)],
        y = ((y + c(-dy2, dy2)) * rt_max)[c(1L, 1L, 2L, 2L)],
        border = NA,
        col = rgb(pal(y), maxColorValue = 255)
      )
    }
    axis(
      side = 4,
      at = par("usr")[3:4],
      labels = c("", ""),
      mgp = c(3, 0.7, 0),
      lwd.ticks = 0
    )
    axis(
      side = 4,
      mgp = c(3, 0.7, 0),
      las = 1,
      lwd = 0,
      lwd.ticks = 1
    )
    text(
      x = par("usr")[2L] + 0 * par("cxy")[1L],
      y = par("usr")[4L] + 0.5 * par("cxy")[2L],
      labels = "growth\nrate,\nper day",
      cex = 0.9,
      adj = c(0, 0),
      xpd = NA
    )
    tdoubling <- c(1:5, 10, 50, 100)
    axis(
      side = 4,
      at = par("usr")[3:4],
      labels = c("", ""),
      mgp = c(3, 4.7, 4),
      lwd.ticks = 0
    )
    axis(
      side = 4,
      at = log(2) / tdoubling,
      labels = tdoubling,
      mgp = c(3, 4.7, 4),
      las = 1,
      lwd = 0,
      lwd.ticks = 1
    )
    text(
      x = par("usr")[2L] + 3.5 * par("cxy")[1L],
      y = par("usr")[4L] + 0.5 * par("cxy")[2L],
      labels = "doubling\ntime,\ndays",
      cex = 0.9,
      adj = c(0, 0),
      xpd = NA
    )

    j <- j + per_plot
  } # loop over plots

  invisible(NULL)
}
