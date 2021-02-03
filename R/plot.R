#' Plot models of epidemic growth
#'
#' @description
#' Methods for plotting `"egf"` objects.
#'
#' @param x
#'   An `"egf"` object returned by [egf()].
#' @param type
#'   A character string indicating a type of plot. Options are
#'   `"interval"` (interval incidence),
#'   `"cumulative"` (cumulative incidence),
#'   `"rt1"`, and `"rt2"` (instantaneous exponential growth rate).
#' @param subset
#'   A named list of atomic vectors with elements specifying levels of
#'   factors in `x$frame`. Only the subset of fitting windows belonging
#'   to these levels is plotted. Use `NULL` (the default) to plot all
#'   fitting windows or if `x$frame` has no factors.
#' @param log
#'   A logical scalar. If `TRUE`, then the dependent variable is plotted
#'   on a logarithmic scale. Unused by `type %in% c("rt1", "rt2")`.
#' @param tol
#'   A non-negative number or `Inf`. Used to define "exceptional"
#'   points when `type = "interval"`. `cases[-1]` is highlighted
#'   according to `control$points_short` if `diff(date) < (1-tol)*m`
#'   and according to `control$points_long` if `diff(date) > (1+tol)*m`,
#'   where `m = median(diff(date))`. In both cases, the value of
#'   `diff(date)` is printed above the point. Assign 0 to highlight
#'   all deviations from `m`. Assign `Inf` to disable highlighting.
#' @param legend
#'   A logical scalar. If `TRUE`, then a legend is displayed in the
#'   right margin. Unused by `type != "rt2"`.
#' @param level
#'   A number in the interval (0,1). The confidence level represented
#'   by confidence intervals on doubling times and confidence bands on
#'   predicted values. Unused by `type = "rt2"`.
#' @param bands
#'   A logical scalar. If `TRUE`, then confidence bands on predicted
#'   values are displayed. Longer run times can be anticipated in this
#'   case. Unused by `type = "rt2"`.
#' @param control
#'   A list of lists defining the appearance of various plot elements,
#'   or otherwise `NULL` (see Details). Unused by `type = "rt2"`.
#' @param per_plot
#'   A positive integer. One plot will display this many time series
#'   when `type = "rt2"`.
#' @param ...
#'   Optional arguments specifying additional graphical parameters
#'   to be recycled for all plots. Currently, only `xlim`, `ylim`,
#'   `ylab`, and `main` are used; see [graphics::plot()]. `xlim`
#'   can be numeric, Date, or character coercible to Date with
#'   `as.Date(xlim)`. If `main` is a character string, then flags
#'   of the form `"%f"`, where `f` is the name of any factor present
#'   in `group_by`, are replaced in a given plot by the factor level
#'   relevant to that plot. Use `"\n"` within `main` to separate
#'   title from subtitle.
#'
#' @return
#' `NULL` (invisibly).
#'
#' @details
#' `control` must be `NULL` or a named list containing some subset
#' of the elements below:
#'
#' \describe{
#' \item{`box`}{
#'   A named list of arguments to [graphics::box()], affecting
#'   the appearance of the box drawn around the plot region.
#'   Only `bty`, `lty`, `lwd`, and `col` are used.
#' }
#' \item{`xax`, `yax`}{
#'   Named lists of arguments to [graphics::axis()], affecting the
#'   appearance of the bottom and left axes. Only `tcl`, `mgp2`,
#'   `col.axis`, `cex.axis`, and `font.axis` are used, where `mgp2`
#'   is the second component of the usual `mgp` argument. The
#'   appearance of the minor (day or month) and major (month or year)
#'   bottom axes can be controlled separately by listing vectors of
#'   length 2. For example, setting
#'   `xax = list(col.axis = c("black", "red"))`
#'   will make the minor axis black and major axis red.
#' }
#' \item{`ylab`, `main`}{
#'   Named lists of arguments to [graphics::title()], affecting
#'   the appearance of axis titles. Only `adj`, `col.lab`, `cex.lab`,
#'   `font.lab`, `col.main`, `cex.main`, and `font.main` are used.
#' }
#' \item{`points`}{
#'   A named list of arguments to [graphics::points()], affecting
#'   the appearance of observed data. Only `pch`, `col`, `bg` and
#'   `cex` are used.
#' }
#' \item{`points_short`, `points_long`}{
#'   Alternatives to `points_main` used to highlight certain points
#'   when `type = "interval"`. See argument `tol`.
#' }
#' \item{`lines`}{
#'   A named list of arguments to [graphics::lines()], affecting
#'   the appearance of predicted curves. Only `lty`, `lwd`, and
#'   `col` are used.
#' }
#' \item{`zero`}{
#'   A named list of arguments to [graphics::lines()], affecting
#'   the appearance of the line drawn at `y = 0` when `type = "rt1"`.
#' }
#' \item{`windows`, `bands`}{
#'   Named lists of arguments to [graphics::polygon()], affecting
#'   the appearance of fitting windows and confidence bands on
#'   predicted incidence curves. Only `col` and `border` are used.
#' }
#' \item{`text_td`, `text_hl`}{
#'   Named lists of arguments to [graphics::text()], affecting
#'   the appearance of text giving initial doubling times and
#'   text above highlighted points. Only `col` and `cex` are
#'   used by `text_td`. `text_hl` also uses `font`, `pos` and
#'   `offset`.
#' }
#' }
#'
#' If `control = NULL`, then it defaults to
#' `get_control_default("plot.egf", type = type)`.
#'
#' If `control` is a list and one of its elements is a partially
#' specified list
#' (e.g., `box = list(bty)` rather than `box = c(bty, lty, lwd, col)`),
#' then that element is filled out with values taken from
#' `get_control_default("plot.egf", type = type)`.
#'
#' If `control` is a list and one if its elements is `NULL`
#' (e.g., `box = NULL`), then the corresponding plot element
#' is suppressed.
#'
#' @export
plot.egf <- function(x,
                     type = c("interval", "cumulative", "rt1", "rt2"),
                     subset = NULL,
                     log = TRUE,
                     tol = 0,
                     legend = FALSE,
                     level = 0.95,
                     bands = FALSE,
                     control = NULL,
                     per_plot = 6L,
                     between_panels = 0.25,
                     ...) {
  type <- match.arg(type)
  interaction0 <- function(...) interaction(..., drop = TRUE, sep = ":")

  group_by <- attr(x$frame, "group_by")
  any_groups <- (length(group_by) > 0L) # TRUE: multiple time series
  any_factors <- (length(x$frame) > 2L) # TRUE: mixed effects model

  ## Augmented frame
  frame_aug <- rbind(x$frame, attr(x$frame, "extra"))
  frame_aug <- cbind(frame_aug, .index = `length<-`(x$index, nrow(frame_aug)))

  ## Reduced frame
  if (any_factors) {
    frame_red <- x$frame[!duplicated(x$index), -(1:2), drop = FALSE]
  }

  ## Subset
  if (any_factors && !is.null(subset)) {
    stop_if_not(
      is.list(subset),
      length(subset) > 0L,
      !is.null(names(subset)),
      m = "`subset` must be a named list or NULL."
    )
    stop_if_not(
      vapply(subset, is.atomic, FALSE),
      lengths(subset) > 0L,
      names(subset) %in% names(frame_red),
      !duplicated(names(subset)),
      unlist(Map(`%in%`, subset, lapply(frame_red[names(subset)], levels))),
      m = "`subset` must specify levels of factors in `x$frame`."
    )
    l <- Map(`%in%`, frame_red[names(subset)], subset)
    w <- Reduce(`&`, l)
    stop_if_not(
      any(w),
      m = "`subset` does not match any fitting windows."
    )
    ## Omit excluded fitting windows
    frame_aug$.index <- factor(frame_aug$.index, levels = levels(x$index)[w])
  }

  ## Split augmented frame by time series
  if (any_groups) {
    ts <- interaction0(frame_aug[group_by])
    frame_aug_split <- split(frame_aug, ts)
  } else {
    frame_aug_split <- list(`1` = frame_aug)
  }

  ## Omit time series now without fitting windows
  omit <- vapply(frame_aug_split, function(d) all(is.na(d$.index)), FALSE)
  frame_aug_split <- frame_aug_split[!omit]

  ## Order time series by date
  frame_aug_split <- lapply(frame_aug_split, function(d) d[order(d[[1L]]), ])

  ## Drop unused levels
  frame_aug_split <- lapply(frame_aug_split, droplevels)

  ## Validate `control` structure (quietly)
  if (is.null(control) || !inherits(control, "list") || is.null(names(control))) {
    control <- get_control_default("plot.egf", type = type)
  } else {
    control_default <- get_control_default("plot.egf", type = type)
    for (pe in intersect(names(control), names(control_default))) {
      if (is.null(control[[pe]])) {
        control_default[pe] <- list(NULL)
      } else if (inherits(control[[pe]], "list") && !is.null(names(control[[pe]]))) {
        s <- intersect(names(control[[pe]]), names(control_default[[pe]]))
        control_default[[pe]][s] <- control[[pe]][s]
      }
    }
    control <- control_default
  }

  if (type %in% c("interval", "cumulative", "rt1")) {
    plot.egf.main(x,
      type = type,
      log = log,
      tol = tol,
      legend = legend,
      bands = bands,
      level = level,
      control = control,
      frame_aug_split = frame_aug_split,
      ...
    )
  } else if (type %in% c("rt2")) {
    plot.egf.heat(x,
      per_plot = per_plot,
      between_panels = between_panels,
      control = control,
      frame_aug_split = frame_aug_split,
      ...
    )
  }
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
  if (type %in% c("interval", "cumulative")) {
    stop_if_not_tf(log)
  } else {
    log <- FALSE
  }
  if (type == "interval") {
    stop_if_not_tf(legend)
  } else {
    legend <- FALSE
  }
  stop_if_not_tf(bands)
  stop_if_not_in_0_1(level)

  ### Set-up ==============================================

  ## A way to avoid conditional `if (type == ...) expr1 else expr2`
  ## in some places
  varname <- switch(type,
    interval = "int_inc",
    cumulative = "cum_inc",
    rt1 = "rt"
  )
  formula <- reformulate("time", varname)

  ## Confidence intervals on doubling times
  show_td <- (
    type %in% c("interval", "rt1") &&
    x$curve %in% c("exponential", "logistic", "richards") &&
    !is.null(control$text_td)
  )
  if (show_td) {
    ci <- confint(x, parm = "tdoubling", level = level, method = "wald")
  }

  ## Axis title (y)
  if (is.null(dots$ylab)) {
    if (type %in% c("interval", "cumulative")) {
      ylab <- sprintf("%s incidence", type)
    } else {
      ylab <- expression("growth rate, day"^{-1})
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
    s <- switch(x$curve, gompertz = "Gompertz", richards = "Richards", x$curve)
    main <- sprintf("Fitted %s model", s)
    if (length(group_by) > 0L) {
      s <- paste(sprintf("%s = %%%s", group_by, group_by), collapse = ", ")
      main <- sprintf("%s\n%s", main, s)
    }
  } else {
    main <- dots$main
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
    if (is.character(main)) {
      ## Replace "%factor_name" with "level_name"
      for (s in group_by) {
        main <- gsub(sprintf("%%%s", s), as.character(d[1L, s]), main, fixed = TRUE)
      }
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
    control$points_main <- control$points
    pty <- c("main", "short", "long")
    if (type == "interval") {
      dt_min <- (1 - tol) * m
      dt_max <- (1 + tol) * m
      dt_enum <- 1L + 1L * (data$dt < dt_min) + 2L * (data$dt > dt_max)
      data$pty <- factor(pty)[dt_enum]
    } else {
      data$pty <- factor("main")
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
    if (!is.null(control$windows)) {
      for (i in seq_len(nlevels(index))) {
        l <- list(
          x = c(t1[i], t2[i], t2[i], t1[i]),
          y = ylim[c(1L, 1L, 2L, 2L)]
        )
        do.call(polygon, c(l, control$windows))
      }
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

    ## Zero line
    if (type == "rt1" && !is.null(control$zero)) {
      do.call(abline, c(list(h = 0), control$zero))
    }

    ## Annotation above exceptional points
    if (type == "interval" && !is.null(control$text_hl)) {
      l <- list(
        formula = formula,
        data = data,
        labels = data$dt,
        subset = (data$pty != pty[1L])
      )
      do.call(text, c(l, control$text_hl))
    }

    ## Confidence band
    if (bands && !is.null(control$bands)) {
      for (p in pred) {
        l <- list(
          x = c(p$time, rev(p$time)),
          y = c(p$lower, rev(p$upper))
        )
        do.call(polygon, c(l, control$bands))
      }
    }

    ## Predicted curve
    if (!is.null(control$lines)) {
      for (p in pred) {
        l <- list(formula = estimate ~ time, data = p)
        do.call(lines, c(l, control$lines))
      }
    }

    ## Doubling time
    if (show_td) {
      y_lu <- add_lines_to_user(0.25, inv_log10(usr[4L], log), log)
      h_lu <- strheight("", cex = control$text_td$cex, font = 1)
      y_e <- add_lines_to_user(0.15, add_height_to_user(h_lu, y_lu, log), log)
      h_e <- strheight("", cex = control$text_td$cex, font = 2)
      for (s in levels(index)) {
        p <- pred[[s]]
        elu <- unlist(ci[match(s, levels(x$index)), length(ci) - 2:0])
        text(
          x = rep.int(mean(p$time[c(1L, nrow(p))]), 2L),
          y = c(y_e, y_lu),
          labels = c(sprintf("%.1f", elu[1L]), sprintf("(%.1f, %.1f)", elu[2L], elu[3L])),
          cex = cex * control$text_td$cex,
          col = control$text_td$col,
          font = c(2, 1),
          adj = c(0.5, 0),
          xpd = NA
        )
      }
      y_lu_again <- add_lines_to_user(0.5, add_height_to_user(h_e, y_e, log), log)
      y_e_again <- add_lines_to_user(0.15, add_height_to_user(h_lu, y_lu_again, log), log)
      text(
        x = rep.int(usr[2L] - 0.5 * strwidth("estimate", cex = control$text_td$cex, font = 2), 2L),
        y = c(y_e_again, y_lu_again),
        labels = c("estimate", "(95% CI)"),
        cex = cex * control$text_td$cex,
        col = control$text_td$col,
        font = c(2, 1),
        adj = c(0.5, 0),
        xpd = NA
      )
      text(
        x = usr[2L],
        y = add_lines_to_user(0.25, add_height_to_user(h_e, y_e_again, log), log),
        labels = "initial doubling time, days:",
        cex = cex * control$text_td$cex,
        col = control$text_td$col,
        font = 1,
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
      if (type %in% c("interval", "cumulative")) {
        l1$at = axTicks(side = 2)
        if (max(l1$at) >= 1e05) {
          l1$labels <- get_labels(l1$at)
          l2$cex.axis <- min(control$yax$cex.axis, get_cex_axis(l1$labels, mex = 0.9 * 4 - control$yax$mgp2))
        }
        do.call(axis, c(l1, l2))
      } else {
        do.call(axis, c(l1, l2))
        if (usr[2L] > 0) {
          axis(
            side = 2,
            line = 4.5,
            at = c(0, usr[4L]),
            labels = c("", ""),
            lwd.ticks = 0
          )
          td <- c(0.5, 1, 5, 10, 50, 100)
          l1$at <- c(log(2) / td)
          l1$labels <- td
          ## FIXME: why 4.5 instead of 4?
          l2$mgp <- c(3, 4 + control$yax$mgp2, 4.5)
          l2$lwd <- 0
          l2$lwd.ticks <- 1
          do.call(axis, c(l1, l2))
        }
      }

      ## Axis title (y)
      if (!is.null(control$ylab)) {
        if (type %in% c("interval", "cumulative")) {
          l <- list(ylab = ylab, line = 4)
          do.call(title, c(l, control$ylab))
        } else {
          l <- list(
            ylab = ylab,
            line = 0.5 + max(strwidth(axTicks(side = 2), units = "inches", cex = control$yax$cex.axis, font = control$yax$font.axis)) / csi + control$yax$mgp2
          )
          do.call(title, c(l, control$ylab))
          if (usr[2L] > 0) {
            rsw <- (strwidth(ylab_again, units = "inches", cex = control$ylab$cex.lab, font = control$ylab$font.lab) * diff(ylim) / pin[2L]) / (0.8 * (usr[2L] - max(0, usr[1L])))
            text(
              ## FIXME: why 7.5 instead of 5?
              x = usr[1L] - cxy[1L] * (7.5 + max(strwidth(td, units = "inches", cex = control$yax$cex.axis, font = control$yax$font.axis)) / csi + control$yax$mgp2),
              y = mean(c(max(0, usr[3L]), usr[4L])),
              labels = ylab_again,
              adj = c(0.5, 0),
              srt = 90,
              xpd = NA,
              col = control$ylab$col.lab,
              cex = control$ylab$cex.lab * par("cex") / max(1, rsw),
              font = control$ylab$font.lab
            )
          }
        }
      }
    }

    ## Axis title (main)
    if (!is.null(control$main)) {
      main_split <- strsplit(main, "\n", fixed = TRUE)[[1L]]
      control$main <- lapply(control$main, rep_len, 2L)
      if (show_td) {
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
      cond <- all(data$dt[data$pty == "main"] == m, na.rm = TRUE)
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
plot.egf.heat <- function(x, per_plot, between_panels,
                          control, frame_aug_split, ...) {
  dots <- list(...)
  stop_if_not_positive_integer(per_plot)

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

  pal <- do.call(colorRamp, control$heat)

  j <- 0L
  while (j < length(frame_aug_split)) {

    layout(matrix(c(seq_len(per_plot), rep.int(per_plot + 1L, per_plot)), ncol = 2L), widths = c(4, 1))
    par(mar = c(0, 1.5, between_panels, 0.5))

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

      text(
        x = par("usr")[1L] + diff(par("usr")[1:2]) * (0.075 * par("pin")[2L] / par("pin")[1L]),
        y = par("usr")[4L] - 0.075 * diff(par("usr")[3:4]),
        labels = names(frame_aug_split)[k],
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

    par(mar = c(0, 1, 0.25, 7))
    plot.new()
    dy2 <- 0.005
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
    td <- c(0.5, 1, 5, 10, 50, 100)
    axis(
      side = 4,
      at = par("usr")[3:4],
      labels = c("", ""),
      mgp = c(3, 4.7, 4),
      lwd.ticks = 0
    )
    axis(
      side = 4,
      at = log(2) / td,
      labels = td,
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
