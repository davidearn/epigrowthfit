#' Plot models of epidemic growth
#'
#' @description
#' Methods for plotting `"egf"` objects.
#'
#' @param x
#'   An `"egf"` object returned by [egf()].
#' @param type
#'   A character string indicating a type of plot. One of
#'   `"interval"` (interval incidence),
#'   `"cumulative"` (cumulative incidence),
#'   `"rt1"`, `"rt2"` (instantaneous exponential growth rate).
#' @param subset
#'   A named list of atomic vectors with elements specifying levels of
#'   factors in `x$frame`. Only the subset of fitting windows belonging
#'   to these levels is plotted. Use `NULL` (the default) to plot all
#'   fitting windows or if `x$frame` has no factors.
#' @param xty
#'   A character string. If `"date"`, then ticks on the time axis are
#'   placed at the start of calendar days, months, or years (depending
#'   on time scale). If `"numeric"`, then time is displayed as a number
#'   of days since an initial date.
#' @param log
#'   A logical scalar. If `TRUE`, then the dependent variable is plotted
#'   on a logarithmic scale. Unused by `type %in% c("rt1", "rt2")`.
#' @param tol
#'   A non-negative number or `Inf`, used to define "exceptional"
#'   points when `type = "interval"`. `cases[-1]` is highlighted
#'   according to `control$points_short` if `diff(date) < (1-tol)*m`
#'   and according to `control$points_long` if `diff(date) > (1+tol)*m`,
#'   where `m = median(diff(date))`. In both cases, the value of
#'   `diff(date)` is printed above the point. Assign 0 to highlight
#'   all deviations from `m`. Assign `Inf` to disable highlighting.
#' @param legend
#'   A logical scalar. If `TRUE`, then a legend is displayed in the
#'   right margin. Unused by `type = "rt2"`.
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
#' @param ...
#'   Optional arguments specifying additional graphical parameters
#'   to be recycled for all plots. Currently, only `xlim`, `ylim`,
#'   `xlab`, `ylab`, and `main` are used; see [graphics::plot()].
#'   `xlim` can be numeric, Date, or character coercible to Date
#'   with `as.Date(xlim)`. If `main` is a character string, then
#'   flags of the form `"%f"`, where `f` is the name of any factor
#'   present in `group_by`, are replaced in a given plot by the
#'   factor level relevant to that plot.
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
#'   Currently, only `bty`, `lty`, `lwd`, and `col` are used.
#' }
#' \item{`xax`, `yax`}{
#'   Named lists of arguments to [graphics::axis()], affecting
#'   the appearance of the bottom and left axes. Currently, only
#'   `tcl`, `mgp2`, `col.axis`, `cex.axis`, and `font.axis` are
#'   used, where `mgp2` is the second component of the usual `mgp`
#'   argument. If `xty = "date"`, then the appearance of the minor
#'   (day or month) and major (month or year) bottom axes can be
#'   controlled separately by listing vectors of length 2.
#'   For example, setting `xax = list(col.axis = c("black", "red"))`
#'   will make the minor axis black and major axis red.
#' }
#' \item{`xlab`, `ylab`, `main`}{
#'   Named lists of arguments to [graphics::title()], affecting
#'   the appearance of axis titles. Currently, only `line`, `adj`,
#'   `col.lab`, `cex.lab`, `font.lab`, `col.main`, `cex.main`,
#'   and `font.main` are used.
#' }
#' \item{`points_main`}{
#'   A named list of arguments to [graphics::points()], affecting
#'   the appearance of observed data. Currently, only `pch`, `col`,
#'   `bg` and `cex` are used.
#' }
#' \item{`points_short`, `points_long`}{
#'   Alternatives to `points_main` used to highlight certain points
#'   when `type = "interval"`. See argument `tol`.
#' }
#' \item{`lines`}{
#'   A named list of arguments to [graphics::lines()], affecting
#'   the appearance of the predicted incidence curve. Currently,
#'   only `lty`, `lwd`, and `col` are used.
#' }
#' \item{`zero`}{
#'   A named list of arguments to [graphics::lines()], affecting
#'   the appearance of the line drawn at `y = 0` when `type = "rt1"`.
#' }
#' \item{`windows`, `bands`}{
#'   Named lists of arguments to [graphics::polygon()], affecting
#'   the appearance of fitting windows and confidence bands on
#'   predicted incidence curves. Currently, only `col` and `border`
#'   are used.
#' }
#' \item{`text_td`, `text_hl`}{
#'   Named lists of arguments to [graphics::text()], affecting
#'   the appearance of text giving doubling times and text above
#'   highlighted points when `type = "interval"`. Currently, only
#'   `adj` (`text_td` only), `pos`, `offset`, `col`, `cex`, and
#'   `font` are used.
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
                     xty = c("date", "numeric"),
                     log = TRUE,
                     tol = 0,
                     legend = FALSE,
                     level = 0.95,
                     bands = FALSE,
                     control = NULL,
                     per_plot = 12L,
                     ...) {
  type <- match.arg(type)
  xty <- match.arg(xty)
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

  if (type %in% c("interval", "cumulative", "rt1")) {
    plot.egf.main(x,
      type = type,
      xty = xty,
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
      xty = xty,
      per_plot = per_plot,
      frame_aug_split = frame_aug_split,
      ...
    )
  }
}

#' @import graphics
#' @importFrom stats reformulate median predict confint
plot.egf.main <- function(x, type, xty, log, tol, legend,
                          bands, level, control, frame_aug_split, ...) {
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
  stop_if_not_tf(legend)
  stop_if_not_tf(bands)
  stop_if_not_in_0_1(level)
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

  ### Set-up ==============================================

  group_by <- attr(x$frame, "group_by")
  any_groups <- (length(group_by) > 0L)

  ## A way to avoid conditional `if (type == ...) expr1 else expr2`
  ## in some places
  varname <- switch(type,
    interval = "int_inc",
    cumulative = "cum_inc",
    rt1 = "rt"
  )
  formula <- reformulate("time", varname)

  ## Confidence intervals on doubling times
  if (type == "interval") {
    ci <- confint(x, parm = "tdoubling", level = level, method = "wald")
  }

  ### Loop over plots =====================================

  op <- par(mar = c(4, 5, 2.7, 1.5 + 5 * legend) + 0.1)
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
    t0 <- 0L
    t1 <- days(d1, since = d0)
    t2 <- days(d2, since = d0)

    data <- data.frame(
      time = days(date, since = d0),
      cum_inc = cumsum(c(0L, cases[-1L])),
      int_inc = cases
    )
    data$dt <- c(NA, diff(data$time))
    data$rt <- c(diff(base::log(cases)) / data$dt[-1L], NA)
    data$rt[!is.finite(data$rt)] <- NA

    ## Predicted curves with confidence bands
    m <- median(data$dt[-1L])
    f <- function(i1, t1, t2, by) {
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
    pred <- Map(f, i1 = i1, t1 = t1, t2 = t2,
                by = switch(type, interval = m, 1))

    ## Axis title (x)
    if (is.null(dots$xlab)) {
      xlab <- switch(xty,
        date = "",
        numeric = sprintf("days since %s", as.character(d0))
      )
    } else {
      xlab <- dots$xlab
    }

    ## Axis title (y)
    if (is.null(dots$ylab)) {
      if (type %in% c("interval", "cumulative")) {
        ylab <- sprintf("%s incidence", type)
      } else {
        ylab <- "instantaneous exponential\ngrowth rate, per day"
      }
    } else {
      ylab <- dots$ylab
    }

    ## Axis title (main)
    if (is.null(dots$main)) {
      s <- switch(x$curve, gompertz = "Gompertz", richards = "Richards", x$curve)
      main <- sprintf("Fitted %s model", s)
      if (any_groups) {
        lv <- as.character(unlist(d[1L, group_by]))
        s <- paste(sprintf("%s = %s", group_by, lv), collapse = ", ")
        main <- sprintf("%s\n(%s)", main, s)
      }
    } else {
      main <- dots$main
      if (is.character(main) && any_groups) {
        for (s in group_by) {
          ## Replace "%factor_name" with "level_name"
          main <- gsub(sprintf("%%%s", s), as.character(d[1L, s]), main, fixed = TRUE)
        }
      }
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
      if (type %in% c("interval", "cumulative")) {
        ymax <- max(data[[varname]], na.rm = TRUE)
        if (log) {
          zero <- ymax^-0.04
          data[[varname]][data[[varname]] == 0L] <- zero
        } else {
          zero <- 0
        }
        ymin <- zero
        ymax <- if (log) ymax^1.04 else ymax * 1.04
        ylim <- c(ymin, ymax)
      } else {
        ylim <- range(data$rt, na.rm = TRUE)
        ylim <- ylim + c(-1, 1) * 0.04 * diff(ylim)
      }
    } else {
      ylim <- dots$ylim
    }

    ## Point styles depending on observation interval
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

    ## Fitting windows
    if (!is.null(control$windows)) {
      for (i in seq_len(nlevels(index))) {
        l <- list(
          x = c(t1[i], t2[i], t2[i], t1[i]),
          y = ylim[c(1, 1, 2, 2)]
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
      l <- list(h = 0)
      do.call(abline, c(l, control$zero))
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
    if (type == "interval" && x$curve %in% c("exponential", "logistic", "richards") && !is.null(control$text_td)) {
      for (s in levels(index)) {
        elu <- unlist(ci[match(s, levels(x$index)), length(ci) - 2:0])
        l <- list(
          labels = sprintf("doubling time:\n%.1f (%.1f, %.1f) days", elu[1L], elu[2L], elu[3L]),
          xpd = NA
        )
        er <- range(pred[[s]]$estimate, na.rm = TRUE)
        if (log) {
          l$y <- er[1L] * 10^(0.25 * diff(log10(er)))
        } else {
          l$y <- er[1L] + 0.25 * diff(er)
        }
        l$x <- min(pred[[s]]$time[pred[[s]]$estimate > l$y], na.rm = TRUE)
        do.call(text, c(l, control$text_td))
      }
    }

    ## Box
    if (!is.null(control$box)) {
      do.call(box, control$box)
    }

    ## Axis (x)
    if (!is.null(control$xax)) {
      if (xty == "date") {
        l <- list(
          left = par("usr")[1L],
          right = par("usr")[2L],
          refdate = d0
        )
        do.call(daxis, c(l, control$xax))
      } else { # "numeric"
        l1 <- list(side = 1L)
        l2 <- lapply(control$xax, "[", 1L)
        l2$mgp <- c(3, l2$mgp2, 0)
        l2$mgp2 <- NULL
        do.call(axis, c(l1, l2))
      }
    }

    ## Axis (y)
    if (!is.null(control$yax)) {
      l1 <- list(
        side = 2L,
        at = axTicks(side = 2L),
        las = 1
      )
      if (max(l1$at) < 1e05) {
        l1$labels <- TRUE
        cex_axis_default <- 0.85
      } else {
        l1$labels <- get_labels(l1$at)
        cex_axis_default <- min(0.85, get_cex_axis(l1$at, mex = 0.9 * control$ylab$line - control$yax$mgp2))
      }
      l2 <- control$yax
      l2$mgp <- c(3, l2$mgp2, 0)
      l2$mgp2 <- NULL
      if (is.na(l2$cex.axis)) {
        l2$cex.axis <- cex_axis_default
      }
      do.call(axis, c(l1, l2))
    }

    ## Axis title (x)
    if (!is.null(control$xlab)) {
      l <- list(xlab = xlab)
      do.call(title, c(l, control$xlab))
    }

    ## Axis title (y)
    if (!is.null(control$ylab)) {
      l <- list(ylab = ylab)
      do.call(title, c(l, control$ylab))
    }

    ## Axis title (main)
    if (!is.null(control$main)) {
      l <- list(main = main)
      do.call(title, c(l, control$main))
    }

    ## Legend
    if (legend) {
      lx <- par("usr")[2L] + 0.02 * diff(par("usr")[1:2])
      ly <- par("usr")[4L] - 0.02 * diff(par("usr")[3:4])
      if (log) {
        ly <- 10^ly
      }
      if (type == "interval") {
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
      } else {
        lexp <- c("obs", NA, NA, "pred")
        lind <- c(TRUE, FALSE, FALSE, TRUE)
      }
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
#' @importFrom grDevices colorRampPalette rgb
plot.egf.heat <- function(x, xty, per_plot, frame_aug_split, ...) {
  dots <- list(...)
  stop_if_not_positive_integer(per_plot)

  group_by <- attr(x$frame, "group_by")
  any_groups <- (length(group_by) > 0L)

  refdate <- do.call(min, lapply(frame_aug_split, `[`, 1L, 1L))
  xlim <- c(0, days(do.call(max, lapply(frame_aug_split, function(d) d[nrow(d), 1L])), since = refdate))
  ylim <- c(0, 1)

  g <- function(d) {
    index <- d$.index
    i12 <- vapply(levels(index), function(s) range(which(index == s)), integer(2L))
    i1 <- i12[1L, ]
    i2 <- i12[2L, ]
    t1 <- days(d[i1, 1L], since = refdate)
    t2 <- days(d[i2, 1L], since = refdate)

    f <- function(i1, t1, t2) {
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
    Map(f, i1 = i1, t1 = t1, t2 = t2)
  }
  pred <- lapply(frame_aug_split, g)
  rt_max <- do.call(max, lapply(do.call(c, pred), function(d) max(d$estimate, na.rm = TRUE)))

  ### Loop over plots =====================================

  op <- par(oma = c(2.5, 0, 3, 0))
  on.exit(par(op))

  pal <- colorRamp(c("#364B9A", "#4A7BB7", "#6EA6CD", "#98CAE1",
                     "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366",
                     "#F67E4B", "#DD3D2D", "#A50026"))

  j <- 0L
  while (j < length(frame_aug_split)) {

    layout(matrix(c(1L, 2L, 3L, 4L, 4L, 4L), nrow = 3L), widths = c(6, 1))
    par(mar = c(0, 3, 0.25, 0.5))

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

      box(lwd = 0.5)
      title(ylab = names(frame_aug_split)[k], line = 1.5)
      if (k == j + 1L) {
        title(
          main = sprintf("Instantaneous exponential growth rate\nper day, by %s", paste(group_by, collapse = ":")),
          line = 0.5,
          adj = 0,
          cex.main = 1.2,
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
      cex.axis = c(0.85, 1.15)
    )

    par(mar = c(0, 1, 0.25, 3))
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
      mgp = c(3, 0.7, 0),
      las = 1,
      col = NA,
      col.ticks = "black"
    )

    j <- j + per_plot
  } # loop over plots

  invisible(NULL)
}
