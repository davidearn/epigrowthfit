#' Plot models of epidemic growth
#'
#' @description
#' Method for plotting `"egf"` objects.
#'
#' @param x
#'   An `"egf"` object returned by [egf()].
#' @param subset
#'   A named list of atomic vectors with elements specifying levels of
#'   factors in `x$frame`. Only the subset of fitting windows belonging
#'   to these levels is plotted. Use `NULL` (the default) to plot all
#'   fitting windows or if `x$frame` has no factors.
#' @param type
#'   A character string indicating the type of incidence plotted.
#' @param log
#'   A logical scalar. If `TRUE`, then incidence is displayed on a
#'   logarithmic scale.
#' @param xty
#'   A character string. If `"date"`, then ticks on the time axis are
#'   placed at the start of calendar days, months, or years (depending
#'   on time scale). If `"numeric"`, then time is displayed as a number
#'   of days since an initial date.
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
#'   right margin.
#' @param bands
#'   A logical scalar. If `TRUE`, then confidence bands on predicted
#'   incidence are displayed. Longer run times should be expected in
#'   this case.
#' @param level
#'   A number in the interval (0,1). The confidence level represented
#'   by confidence intervals on doubling times and confidence bands on
#'   predicted incidence.
#' @param control
#'   A list of lists defining the appearance of various
#'   plot elements, or otherwise `NULL` (see Details).
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
#'   controlled separately by listing vectors of length 2. For example,
#'   setting `xax = list(col.axis = c("black", "red"))` will make
#'   the minor axis black and major axis red.
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
#' \item{`windows`, `bands`}{
#'   Named lists of arguments to [graphics::polygon()], affecting
#'   the appearance of fitting windows and confidence bands on
#'   predicted incidence curves. Currently, only `col` and `border`
#'   are used.
#' }
#' \item{`text_td`, `text_hl`}{
#'   Named lists of arguments to [graphics::text()], affecting
#'   the appearance of text giving doubling times and text above
#'   highlighted points. Currently, only `adj` (`text_td` only),
#'   `pos`, `offset`, `col`, `cex`, and `font` are used.
#' }
#' }
#'
#' If `control = NULL`, then it defaults to
#' `get_control_default("plot.egf")`.
#'
#' If `control` is a list and one of its elements is a partially
#' specified list
#' (e.g., `box = list(bty)` rather than `box = c(bty, lty, lwd, col)`),
#' then that element is filled out with values taken from
#' `get_control_default("plot.egf")`.
#'
#' If `control` is a list and one if its elements is `NULL`
#' (e.g., `box = NULL`), then the corresponding plot element
#' is suppressed.
#'
#' @export
#' @import graphics
#' @importFrom stats reformulate median predict confint
plot.egf <- function(x,
                     subset = NULL,
                     type = c("interval", "cumulative"),
                     log = TRUE,
                     xty = c("date", "numeric"),
                     tol = 0,
                     legend = FALSE,
                     bands = FALSE,
                     level = 0.95,
                     control = NULL,
                     ...) {

  ### Argument validation #################################

  any_factors <- (length(x$frame) > 2L)
  if (any_factors) {
    ## Reduced frame
    frame_red <- x$frame[!duplicated(x$index), -(1:2), drop = FALSE]
    group_by <- attr(x$frame, "group_by")
    any_groups <- (length(group_by) > 0L)

    if (is.null(subset)) {
      index_levels_for_plot <- levels(x$index)
    } else {
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
        unlist(Map("%in%", subset, lapply(frame_red[names(subset)], levels))),
        m = "`subset` must specify levels of factors in `x$frame`."
      )
      w <- Reduce("&", Map("%in%", frame_red[names(subset)], subset))
      stop_if_not(
        any(w),
        m = "`subset` does not match any fitting windows."
      )
      index_levels_for_plot <- levels(x$index)[w]

      if (any_groups) {
        i <- which(names(subset) %in% group_by)
        w <- Reduce("&", Map("%in%", frame_red[names(subset)[i]], subset[i]))
        frame_red <- frame_red[w, , drop = FALSE]
      }
    }
  } else {
    any_groups <- FALSE
    index_levels_for_plot <- levels(x$index)
  }

  type <- match.arg(type)
  stop_if_not_tf(log)
  xty <- match.arg(xty)
  stop_if_not_tf(legend)
  if (type == "interval") {
    stop_if_not(
      is.numeric(tol),
      length(tol) == 1L,
      tol >= 0,
      m = "`tol` must be a non-negative number or `Inf`."
    )
  }
  stop_if_not_tf(bands)
  stop_if_not_in_0_1(level)

  if (is.null(control) || !inherits(control, "list") || is.null(names(control))) {
    control <- get_control_default("plot.egf")
  } else {
    control_default <- get_control_default("plot.egf")
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
  dots <- list(...)

  ### Setting up loop over plots ##########################

  ## Augmented frame
  frame_aug <- rbind(x$frame, attr(x$frame, "extra"))
  frame_aug <- cbind(
    index = `length<-`(x$index, nrow(frame_aug)),
    frame_aug
  )

  interaction0 <- function(...) interaction(..., drop = TRUE, sep = ":")

  if (any_factors) {
    ## Split the augmented frame by interaction
    frame_aug_split <- split(frame_aug, interaction0(frame_aug[-(1:3)]))
    ## Split the reduced frame by group of interactions
    if (length(group_by) > 0L) {
      frame_red_split <- split(frame_red, interaction0(frame_red[group_by]))
    } else {
      frame_red_split <- list("1" = frame_red)
    }
    ## Merge grouped interactions in split augmented frame
    frame_aug_split <- lapply(frame_red_split, function(d) {
      do.call(rbind, frame_aug_split[as.character(interaction0(d))])
    })
  } else {
    frame_aug_split <- list("1" = frame_aug)
  }
  ## Order by date
  frame_aug_split <- lapply(frame_aug_split, function(d) d[order(d[[2L]]), ])

  stop_if_not(
    vapply(frame_aug_split, function(d) all(diff(d[[2L]]) > 0), FALSE),
    m = paste0(
      "Plots of multiple time series in one panel are\n",
      "not supported (yet). See `group_by` details in\n",
      "`help(\"plot.egf\")`."
    )
  )
  if (type == "cumulative") {
    stop_if_not(
      vapply(frame_aug_split, function(d) !anyNA(d[-1L, 3L]), FALSE),
      m = paste0(
        "Missing values in interval incidence time series\n",
        "preventing calculation of cumulative incidence."
      )
    )
  }

  ## A way to avoid conditional `if (type == ...) expr1 else expr2`
  ## in some places
  varname <- sprintf("%s_inc", substr(type, start = 1L, stop = 3L))
  formula <- reformulate("time", varname)

  ## Confidence intervals on doubling times
  ci <- confint(x, parm = "tdoubling", level = level, method = "wald")

  ### Loop over plots #####################################

  op <- par(
    mar = c(4, 5, 2.7, 0.5 + 6 * legend) + 0.1,
    bty = "l",
    xaxs = "i",
    yaxs = "i",
    las = 1
  )
  on.exit({
    assign("egf.par", par("mar", "plt"), envir = .epigrowthfit)
    par(op)
  })

  for (d in frame_aug_split) {

    ### Setting up plot ===================================

    index <- droplevels(factor(d[[1L]], levels = index_levels_for_plot))
    date <- d[[2L]]
    cases <- d[[3L]]
    i12 <- vapply(levels(index), function(s) range(which(index == s)), integer(2L))
    d0 <- date[1L]
    d1 <- date[i12[1L, ]]
    d2 <- date[i12[2L, ]]
    t0 <- 0L
    t1 <- days(d1, since = d0)
    t2 <- days(d2, since = d0)

    data <- data.frame(
      time = days(date, since = d0),
      cum_inc = cumsum(c(0L, cases[-1L])),
      int_inc = cases,
      dt = c(NA, ddiff(date))
    )

    ## A way to artificially include zeros on logarithmic scale
    ymax <- max(data[[varname]], na.rm = TRUE)
    if (log) {
      zero <- ymax^-0.04
      data[[varname]][data[[varname]] == 0L] <- zero
    } else {
      zero <- 0
    }

    ## Predicted curves with confidence bands
    m <- median(data$dt, na.rm = TRUE)
    f <- function(i1, t1, t2, by) {
      p <- predict(x,
        subset = if (length(d) > 3L) d[i1, -(1:3), drop = FALSE],
        time = seq.int(from = 0L, to = t2 - t1, by = by),
        se = bands
      )
      if (bands) {
        p <- confint(p, level = level, log = FALSE)
      }
      p <- p[[varname]]
      p[[1L]] <- t1 + p[[1L]]
      if (type == "cumulative" && i1 > 1L) {
        p[-1L] <- sum(cases[2L:i1]) + p[-1L]
      }
      p
    }
    pred <- Map(f, i1 = i12[1L, ], t1 = t1, t2 = t2,
                by = if (type == "cumulative") 1L else m)

    ## Axis title (x)
    if (is.null(dots$xlab)) {
      xlab <- switch(xty,
        date = "date",
        numeric = sprintf("days since %s", as.character(d0))
      )
    } else {
      xlab <- dots$xlab
    }

    ## Axis title (y)
    if (is.null(dots$ylab)) {
      ylab <- sprintf("%s incidence", type)
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
      ymin <- zero
      ymax <- if (log) ymax^1.04 else ymax * 1.04
      ylim <- c(ymin, ymax)
    } else {
      ylim <- dots$ylim
    }

    ## Point styles depending on observation interval
    dt_min <- (1 - tol) * m
    dt_max <- (1 + tol) * m
    dt_enum <- 1L + (type == "interval") *
      (1L * (data$dt < dt_min) + 2L * (data$dt > dt_max))
    pty <- c("main", "short", "long")
    data$pty <- factor(pty[dt_enum], exclude = NA)

    ### Plotting ==========================================

    plot.new()
    plot.window(xlim = xlim, ylim = ylim, log = if (log) "y" else "")

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
        at = axTicks(side = 2L)
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
    if (!is.null(control$text_td) && type == "interval" && x$curve %in% c("exponential", "logistic", "richards")) {
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
      if (type == "cumulative") {
        lexp <- c("obs", NA, NA, "pred")
        lind <- c(TRUE, FALSE, FALSE, TRUE)
      } else { # interval
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
      }
      get_el <- function(pe) sapply(control[sprintf("points_%s", pty)], "[[", pe)
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
        lty = c(rep(NA, 3L), control$lines$lty)[lind],
        lwd = c(rep(NA, 3L), control$lines$lwd)[lind],
        col = c(get_el("col"), control$lines$col)[lind]
      )
    }
  }

  invisible(NULL)
}
