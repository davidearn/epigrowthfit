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
#'   `type != "rt2"` displays one time series per plot.
#'   `type = "rt2"` displays `per_plot` time series per plot.
#' @param subset
#'   An expression to be evaluated in the combined model frame.
#'   Must evaluate to a logical vector or list of logical vectors
#'   indexing rows of the model frame, and thus fitting windows.
#'   Only indexed fitting windows are plotted. The default (`NULL`)
#'   is to plot all fitting windows.
#' @param order
#'   An expression to be evaluated in the combined model frame,
#'   typically a call to [order()], determining the order in which
#'   time series are plotted. The default (`NULL`) is equivalent
#'   to `seq_len(nrow(frame))`.
#' @param cache
#'   An `"egf_plot_cache"` object returned by a previous evaluation
#'   of `plot.egf(x)`. Predicted values and standard errors stored
#'   in `cache` will be reused (rather than recomputed) if relevant
#'   to the current function call (i.e., to `type` and `subset`).
#' @param plot
#'   A logical scalar. If `FALSE`, then nothing is plotted.
#'   Useful when only the return value is desired.
#' @param time_as
#'   A character string indicating how time is displayed on the
#'   bottom axis. The options are: as a calendar (`"Date"`) and
#'   as a number of days since the earliest time point (`"numeric"`).
#' @param log
#'   A logical scalar. If `TRUE`, then the dependent variable is
#'   plotted on a logarithmic scale. Unused by `type = "rt[12]"`.
#' @param show_fits
#'   A logical scalar. If `TRUE`, then predicted curves are drawn.
#'   Unused by `type = "rt2"`.
#' @param show_bands
#'   A logical scalar. If `TRUE`, then confidence bands on
#'   predicted curves are drawn. Unused by `type = "rt2"`.
#' @param show_tdoubling
#'   A logical scalar. If `TRUE`, then initial doubling time
#'   estimates and corresponding confidence intervals are printed
#'   in the top margin. Unused by `type = "rt2"`. Unsupported
#'   for `x$curve` not in `c("exponential", "logistic", "richards")`.
#' @param show_legend
#'   A logical scalar. If `TRUE`, then a legend is displayed in
#'   the right margin. Unsupported for `type != "interval"`.
#' @param level
#'   A number in the interval (0,1). The confidence level desired
#'   when `show_tdoubling = TRUE` or `show_bands = TRUE`.
#' @param per_plot
#'   A positive integer. The number of time series displayed
#'   in one plot. Unsupported for `type != "rt2"`.
#' @param control
#'   A named list controlling the appearance of various plot elements
#'   (see Details). The default (`NULL`) is to use
#'   `get_control_default("plot.egf", type = type, time_as = time_as)`.
#' @param xlim,ylim
#'   Optional numeric vectors specifying axis limits, which are
#'   recycled for all plots. If `time_as = "Date"`, then `xlim`
#'   can instead be a Date vector or a character vector coercible
#'   to Date via `as.Date(xlim)`. `ylim` is unused by `type = "rt2"`.
#' @param main,sub,xlab,ylab,ylab_outer,plab
#'   Optional character strings or expressions used to generate plot
#'   (`main`, `sub`), axis (`xlab`, `ylab`, `ylab_outer`), and panel
#'   (`plab`) titles. `sub` is unused by `type = "rt2"`. `ylab_outer`
#'   is used only by `type = "rt[12]"`. `plab` is used only by
#'   `type = "rt2"`. When `type != "rt2"`, `main` and `sub` are
#'   evaluated in the combined model frame to generate unique
#'   (sub)titles for each plot. When `type = "rt2"`, `plab` is
#'   evaluated similarly to generate unique titles for each panel.
#'   Note that [`plotmath`][grDevices::plotmath] expressions are
#'   not supported for `main`, `sub`, and `plab` in these cases.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' An `"egf_plot_cache"` object. If `cache` was used, then this object
#' is the result of augmenting `cache` with any new computations.
#'
#' Note that this object is returned _in spite of_ any errors thrown
#' during plotting, to avoid wasting computation time.
#'
#' @details
#' The combined model frame is `do.call(cbind, unname(x$frame_par))`.
#' If a variable occurs more than once there, because it appears
#' in multiple model frames, then only the earliest instance is
#' retained. Except in unusual cases, all instances of a variable
#' are identical, and no information is lost.
#'
#' `plot.egf()` will _not_ detect mismatch between `x` and `cache`.
#' In other words, constructions like `plot(x2, cache = plot(x1))`
#' should _not_ be expected to produce correct results.
#'
#' `control` is a named list defining some subset of the available
#' control parameters and their respective components (see section
#' "Control parameters"). A list of modifiable options for each
#' control parameter and their default values can be obtained with
#' `get_control_default("plot.egf", type = type, time_as = time_as)`.
#' Unspecified options take their values in this list. Unsupported
#' options are silently discarded. To suppress a plot element, set
#' the corresponding control parameter to `NULL`,
#' as in `control = list(box = NULL)`.
#'
#' # Control parameters
#'
#' For `type != "rt2"`:
#' \describe{
#' \item{`box`}{
#'   A named list of arguments to [graphics::box()],
#'   affecting the appearance of the box drawn around the plot region.
#' }
#' \item{`axis`}{
#'   A named list of the form `list(x, y)` affecting the appearance of
#'   the bottom and left axes, respectively.
#'   If `time_as = "Date"`, in which case there are minor (day or month)
#'   and major (month or year) bottom axes, then `x` should be a named
#'   list of the form `list(minor, major)`, and `minor`, `major`, and
#'   `y` should all be named lists of arguments to [graphics::axis()].
#'   If `time_as = "numeric"`, in which case there is one bottom axis,
#'   then `x` should itself be a named list of arguments to
#'   [graphics::axis()].
#' }
#' \item{`title`}{
#'   A named list of the form `list(main, sub, xlab, ylab)`, affecting
#'   the appearance of plot and axis titles. `main`, `sub`, `xlab`, and
#'   `ylab` should all be lists of arguments to [graphics::title()].
#'   Note that subtitles are handled specially: they are placed under
#'   main titles in the top margin, rather than in the bottom margin
#'   (the usual behavior of [graphics::title()]).
#' }
#' \item{`rect`}{
#'   A named list of arguments to [graphics::rect()],
#'   affecting the appearance of fitting windows.
#' }
#' \item{`polygon`}{
#'   A named list of arguments to [graphics::polygon()],
#'   affecting the appearance of confidence bands on predicted curves.
#' }
#' \item{`lines`}{
#'   A named list of arguments to [graphics::lines()],
#'   affecting the appearance of predicted curves.
#' }
#' \item{`points`}{
#'   A named list of arguments to [graphics::points()],
#'   affecting the appearance of observed data.
#' }
#' \item{`tdoubling`}{
#'   A named list of the form `list(caption, estimate, ci)`,
#'   affecting the appearance of initial doubling times printed
#'   when `x$curve %in% c("exponential", "logistic", "richards")`.
#'   `caption`, `estimate`, and `ci` should each be named lists
#'   of arguments to [graphics::text()].
#' }
#' }
#' For `type = "interval"`:
#' \describe{
#' \item{`special`}{
#'   A named list of the form
#'   `list(tol, points = list(short, long), text)`,
#'   affecting the appearance of "exceptional" points.
#'   `tol` should be a non-negative number or `Inf`. `short` and `long`
#'   should be named lists of arguments to [graphics::points()], and
#'   `text` should be a named list of arguments to [graphics::text()].
#'   Within a time series `data.frame(time, x)`, `x[-1]` is highlighted
#'   according to `short` if `diff(time) < (1-tol)*m` and
#'   according to `long`  if `diff(time) > (1+tol)*m`,
#'   where `m = median(diff(time))`.
#'   In both cases, the value of `diff(time)` is printed above the point
#'   according to `text`. Use `tol = 0` to highlight all deviations from
#'   `m`. Use `tol = Inf` to disable highlighting.
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
#'   A named list of arguments to [grDevices::colorRamp()]
#'   defining the heat map's color palette.
#' }
#' \item{`ips`}{
#'   A non-negative number defining the space between panels
#'   as a number of margin lines.
#' }
#' }
#'
#' @export
#' @importFrom stats median fitted predict complete.cases
plot.egf <- function(x,
                     type = c("interval", "cumulative", "rt1", "rt2"),
                     subset = NULL,
                     order = NULL,
                     cache = NULL,
                     plot = TRUE,
                     time_as = c("Date", "numeric"),
                     log = TRUE,
                     show_fits = TRUE,
                     show_bands = FALSE,
                     show_tdoubling = FALSE,
                     show_legend = FALSE,
                     level = 0.95,
                     per_plot = 6L,
                     control = NULL,
                     xlim = NULL,
                     ylim = NULL,
                     main = NULL,
                     sub = NULL,
                     xlab = NULL,
                     ylab = NULL,
                     ylab_outer = NULL,
                     plab = NULL,
                     ...) {
  ## FIXME: `order` _actually_ permutes `seq_along(levels(window))`
  ## instead of `seq_along(levels(ts))`, leading to indexing gore ...

  type <- match.arg(type)
  stop_if_not_true_false(plot)
  if (type == "rt2") {
    show_fits <- TRUE
    show_bands <- FALSE
    show_tdoubling <- FALSE
  } else {
    stop_if_not_true_false(show_fits)
    stop_if_not_true_false(show_bands)
    stop_if_not_true_false(show_tdoubling)
    show_tdoubling <- show_tdoubling &&
      x$curve %in% c("exponential", "logistic", "richards")
    if (show_tdoubling || show_bands) {
      stop_if_not_number_in_interval(level, 0, 1, "()")
    }
  }
  frame_par <- do.call(cbind, unname(x$frame_par))
  frame_par <- frame_par[!duplicated(names(frame_par))]
  subset <- subset_to_index(substitute(subset), frame_par, parent.frame())
  stop_if_not(
    length(subset) > 0L,
    m = "`subset` must index at least one fitting window."
  )

  ## This code _could_ be run _after_ augmenting `cache`,
  ## but running it _before_ avoids waste of computation
  ## time in the event of errors. This really only matters
  ## for users who have not assigned `plot(x)`.
  if (plot) {
    order <- order_to_index(substitute(order), frame_par, parent.frame())
    if (type == "rt2") {
      plab <- label_to_character(substitute(plab), frame_par, parent.frame())
    } else {
      main <- label_to_character(substitute(main), frame_par, parent.frame())
      sub <- label_to_character(substitute(sub), frame_par, parent.frame())
    }

    if (!is.null(xlim)) {
      if (is.character(xlim)) {
        xlim <- as.Date(xlim)
      }
      stop_if_not(
        is.numeric(xlim) || inherits(xlim, "Date"),
        length(xlim) == 2L,
        xlim[1L] < xlim[2L],
        m = "Invalid `xlim`."
      )
    }
    if (!is.null(ylim)) {
      stop_if_not(
        is.numeric(ylim),
        length(ylim) == 2L,
        ylim[1L] < ylim[2L],
        m = "Invalid `ylim`."
      )
    }

    time_as <- match.arg(time_as)
    if (type == "rt2") {
      stop_if_not_positive_integer(per_plot)
    } else {
      if (type == "rt1") {
        log <- FALSE
      } else {
        stop_if_not_true_false(log)
      }
      if (type == "interval") {
        stop_if_not_true_false(show_legend)
      } else {
        show_legend <- FALSE
      }
    }

    frame_ts <- x$frame_ts
    wl <- as.character(x$endpoints$window)
    tsl <- as.character(x$endpoints$ts)

    wl_subset <- wl[subset]
    frame_ts$window <- factor(frame_ts$window, levels = wl_subset)

    tsl_subset <- intersect(tsl[order], frame_ts$ts[!is.na(frame_ts$window)])
    frame_ts$ts <- factor(frame_ts$ts, levels = tsl_subset)

    if (type == "cumulative") {
      keep <- !tapply(frame_ts$x, frame_ts$ts, function(x) anyNA(x[-1L]))
      if (any(!keep)) {
        warning(
          "Missing values preventing calculation of\n",
          "cumulative incidence. These time series\n",
          "will not be displayed:\n\n",
          paste0("  ", tsl_subset[!keep], collapse = "\n"),
          call. = FALSE
        )
        tsl_subset <- tsl_subset[keep]
        frame_ts$ts <- factor(frame_ts$ts, levels = tsl_subset)
        subset <- subset[tsl[subset] %in% tsl_subset]
        wl_subset <- wl[subset]
        frame_ts$window <- factor(frame_ts$window, levels = wl_subset)
      }
      if (all(!keep)) {
        stop("Nothing left to do ...", call. = FALSE)
      }
    }
  }

  ## Initialize `cache` if not supplied
  if (is.null(cache)) {
    cache <- data.frame(
      var = factor(),
      ts = factor(),
      window = factor(),
      time = numeric(0L),
      estimate = numeric(0L),
      se = numeric(0L)
    )
  } else {
    stop_if_not(
      inherits(cache, "egf_plot_cache"),
      m = "`cache` must inherit from class \"egf_plot_cache\"."
    )
  }

  ## If necessary, augment with fitted values of "log_r" and
  ## standard errors. Since `x` already stores both of these,
  ## no computation is required and it makes sense to include
  ## everything regardless of `show_tdoubling` and `subset`.
  cn <- names(cache)
  if (x$curve %in% c("exponential", "logistic", "richards") &&
      !"log_r" %in% levels(cache$var)) {
    ft <- fitted(x, par = "log_r", link = TRUE, se = TRUE)
    names(ft)[names(ft) == "par"] <- "var"
    ft$time <- NA_real_
    cache <- rbind(cache, ft[cn])
  }

  ## If necessary, augment with predicted values of whatever
  ## is being plotted and standard errors. This _could_ require
  ## a lot of computation, so we are careful here to do exactly
  ## what is necessary.
  if (show_fits || show_bands) {
    what <- switch(type, interval = "log_int_inc", cumulative = "log_cum_inc", "log_rt")
    window <- x$frame_ts$window
    ok <- cache$var == what & !(show_bands & is.na(cache$se))
    need <- setdiff(levels(window)[subset], cache$window[ok])
    if (length(need) > 0L) {
      if (type == "interval") {
        time <- x$frame_ts$time
        dt <- tapply(time, factor(window, levels = need), function(x) median(diff(x)))
      } else {
        dt <- 1
      }
      k <- match(need, levels(window), 0L)
      time_split <- Map(seq.int,
        from = x$endpoints$start[k],
        to = x$endpoints$end[k],
        by = dt
      )
      pd <- predict(x,
        what = what,
        time = unlist(time_split),
        window = rep.int(x$endpoints$window[k], lengths(time_split)),
        log = TRUE,
        se = show_bands
      )
      if (!show_bands) {
        pd$se <- NA_real_
      }
      cache <- rbind(cache, pd[cn])
      ord <- do.call(base::order, cache[cn[c(1L, 3:4, 6L)]])
      cache <- cache[ord, , drop = FALSE]
      keep <- !duplicated(cache[cn[c(1L, 3:4)]])
      cache <- cache[keep, , drop = FALSE]
    }
  }
  class(cache) <- c("egf_plot_cache", "data.frame")

  ## If not plotting, then return `cache`
  if (!plot) {
    return(invisible(cache))
  }

  ## If plotting, then create an instruction to return
  ## `cache` if the low level plot function throws an error
  on.exit({
    message("Augmented `cache` returned despite error ...")
    return(invisible(cache))
  })

  cache1 <- cache
  if (nrow(cache1) > 0L) {
    ## Extract only those parts of `cache` needed by the
    ## low level plot function
    cache1[cn[1:3]] <- Map(factor,
      x = cache1[cn[1:3]],
      levels = list(c(if (type != "rt2") "log_r", what), tsl_subset, wl_subset)
    )
    keep <- complete.cases(cache1[cn[1:3]])
    cache1 <- cache1[keep, , drop = FALSE]
    ord <- do.call(base::order, cache1[cn[1:3]])
    cache1 <- cache1[ord, , drop = FALSE]

    ## Compute confidence intervals as needed
    if (type != "rt2"  && (show_tdoubling || show_bands)) {
      l <- (show_tdoubling & cache1$var == "log_r") | (show_bands & cache1$var == what)
      cache1[l, c("lower", "upper")] <- do_wald(
        estimate = cache1$estimate[l],
        se = cache1$se[l],
        level = level
      )
    }
  }

  ## Last few arguments to low level plot function
  origin <- attr(x$frame_ts, "origin")
  control_default <- get_control_default("plot.egf", type = type, time_as = time_as)
  if (is.null(control)) {
    control <- control_default
  } else {
    control <- clean(control, template = control_default)
  }

  ## Need to subset and order plot/axis/panel titles so that they
  ## correspond elementwise to `tsl_subset`
  subset <- order[order %in% subset]
  m <- match(tsl_subset, tsl[subset], 0L)

  if (type == "rt2") {
    do_heat_plot(
      cache = cache1,
      origin = origin,
      time_as = time_as,
      per_plot = per_plot,
      control = control,
      xlim = xlim,
      main = main,
      xlab = xlab,
      ylab = ylab,
      ylab_outer = ylab_outer,
      plab = plab[subset][m]
    )
  } else {
    do_curve_plot(
      frame_ts = frame_ts,
      cache = cache1,
      origin = origin,
      curve = x$curve,
      type = type,
      time_as = time_as,
      log = log,
      show_fits = show_fits,
      show_bands = show_bands,
      show_tdoubling = show_tdoubling,
      show_legend = show_legend,
      level = level,
      control = control,
      xlim = xlim,
      ylim = ylim,
      main = main[subset][m],
      sub = sub[subset][m],
      xlab = xlab,
      ylab = ylab,
      ylab_outer = ylab_outer
    )
  }

  ## Discard exit instructions when low level plot function
  ## runs without error
  on.exit()
  invisible(cache)
}

#' @import graphics
#' @importFrom stats median
do_curve_plot <- function(frame_ts, cache,
                          origin, curve, type, time_as, log,
                          show_fits, show_bands, show_tdoubling, show_legend,
                          level, control, xlim, ylim,
                          main, sub, xlab, ylab, ylab_outer) {
  ### Set up ==============================================

  N <- nlevels(frame_ts$ts)
  frame_ts_split <- split(frame_ts, frame_ts$ts)
  cache_split <- split(cache, cache$ts)
  inv_log <- if (log) function(x) 10^x else identity
  formula <- as.formula(call("~", as.name(type), as.name("time")))
  what <- switch(type, interval = "log_int_inc", cumulative = "log_cum_inc", "log_rt")

  xlim_bak <- xlim
  ylim_bak <- ylim
  xlab_bak <- xlab

  ## Plot (sub)title
  if (is.null(main)) {
    if (curve %in% c("gompertz", "richards")) {
      substr(curve, 1L, 1L) <- toupper(substr(curve, 1L, 1L))
    }
    main <- rep_len(sprintf("Fitted %s model", curve), N)
    sub <- names(frame_ts_split)
  }

  ## Axis title (y)
  if (is.null(ylab)) {
    if (type == "rt1") {
      ylab <- "per capita growth rate, per day"
    } else {
      ylab <- paste(type, "incidence")
    }
  }
  if (is.null(ylab_outer)) {
    ylab_outer <- "doubling time, days"
  }

  ## Plot margins
  mar <- switch(type,
    rt1 = c(3, 4 + (is.null(ylim) || ylim[2L] > 0) * 4, 4, 1) + 0.1,
    c(3, 5, 4, 1 + 5.5 * show_legend) + 0.1
  )
  op <- par(mar = mar)
  on.exit(par(op))


  ### Loop over plots =====================================

  for (k in seq_len(N)) {
    ### Set up for plot ===================================

    frame_ts <- droplevels(frame_ts_split[[k]])
    shift <- min(frame_ts$time)
    wl <- levels(frame_ts$window)

    cache <- droplevels(cache_split[[k]])
    cache_log_r <- cache[cache$var == "log_r", , drop = FALSE]
    cache_what <- cache[cache$var == what, , drop = FALSE]
    cache_what$time <- cache_what$time - shift
    cache_what_split <- split(cache_what, cache_what$window)

    data <- data.frame(
      time = frame_ts$time - shift,
      dt = c(NA, diff(frame_ts$time))
    )
    data[[type]] <- switch(type,
      interval = c(NA, frame_ts$x[-1L]),
      cumulative = cumsum(c(0, frame_ts$x[-1L])),
      rt1 = c(NA, diff(base::log(frame_ts$x))) / data$dt
    )
    data[[type]][!is.finite(data[[type]])] <- NA

    t12 <- tapply(data$time, frame_ts$window, range, simplify = FALSE)
    t1 <- vapply(t12, `[`, 0, 1L)
    t2 <- vapply(t12, `[`, 0, 2L)

    ## Axis limits (x)
    if (is.null(xlim_bak)) {
      xlim <- c(0, max(data$time) * 1.04)
    } else {
      if (inherits(xlim_bak, "Date")) {
        xlim <- julian(xlim_bak, origin = origin + shift)
      }
    }

    ## Axis limits (y)
    if (is.null(ylim_bak)) {
      if (type == "rt1") {
        ylim <- range(data[[type]], na.rm = TRUE)
        ylim[1L] <- max(-base::log(2), ylim[1L])
        ylim[2L] <- min( base::log(2), ylim[2L])
        ylim <- ylim + c(-1, 1) * 0.04 * (ylim[2L] - ylim[1L])
      } else {
        ymax <- max(data[[type]], na.rm = TRUE)
        if (log) {
          zero <- ymax^-0.04
          data[[type]][data[[type]] == 0] <- zero
        } else {
          zero <- 0
        }
        ylim <- c(zero, if (log) ymax^1.04 else ymax * 1.04)
      }
    }

    ## Axis title (x)
    if (is.null(xlab_bak)) {
      xlab <- switch(time_as,
        Date = "",
        sprintf("number of days since %s", origin + shift)
      )
    }

    ## Point highlighting according to observation interval
    ## (relative to median)
    m <- median(data$dt[-1L])
    pty <- c("median", "short", "long")
    if (type == "interval" && !is.null(control$special)) {
      tol <- control$special$tol
      dt_min <- (1 - tol) * m
      dt_max <- (1 + tol) * m
      dt_enum <- 1L + 1L * (data$dt < dt_min) + 2L * (data$dt > dt_max)
      data$pty <- factor(pty)[dt_enum]
    } else {
      data$pty <- factor("median")
    }
    control$special$points$median <- control$points # hack

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
        ybottom = inv_log(usr[3L]),
        xright = t2,
        ytop = inv_log(usr[4L])
      )
      do.call(rect, c(l, control$rect))
    }

    ## Observed data
    for (s in pty) {
      if (!is.null(control$special$points[[s]])) {
        l <- list(
          formula = formula,
          data = data,
          subset = (data$pty == s)
        )
        do.call(points, c(l, control$special$points[[s]]))
      }
    }

    ## Annotation above exceptional points
    if (type == "interval" && !is.null(control$special$text)) {
      l <- list(
        formula = formula,
        data = data,
        labels = data$dt,
        subset = (data$pty != pty[1L])
      )
      do.call(text, c(l, control$special$text))
    }

    ## Confidence bands
    if (show_bands && !is.null(control$polygon)) {
      for (pd in cache_what_split) {
        l <- list(
          x = c(pd$time, rev(pd$time)),
          y = c(exp(pd$lower), rev(exp(pd$upper)))
        )
        do.call(polygon, c(l, control$polygon))
      }
    }

    ## Predicted curves
    if (show_fits && !is.null(control$lines)) {
      for (pd in cache_what_split) {
        pd$estimate <- exp(pd$estimate)
        if (type == "cumulative") {
          c0 <- data[[type]][match(pd$time[1L], data$time, 0L)]
          pd$estimate <- c(0, c0 + pd$estimate[-1L])
        }
        l <- list(formula = estimate ~ time, data = pd)
        do.call(lines, c(l, control$lines))
      }
    }

    ## Asymptotes
    if (type == "rt1") {
      if (!is.null(control$abline)) {
        l <- list(h = 0)
        do.call(abline, c(l, control$abline))
      }
      if (!is.null(control$segments)) {
        r <- exp(cache_log_r$estimate)
        l <- list(x0 = t1, y0 = r, x1 = t2, y1 = r)
        do.call(segments, c(l, control$segments))
      }
    }

    ## Initial doubling times
    if (show_tdoubling && !is.null(control$tdoubling)) {
      ## Much ado about finding the right user coordinates
      ## when log scale is in effect
      ciy <- add_lines_to_user(0.25, inv_log(usr[4L]), log)
      cih <- strheight("",
        units = "user",
        cex = control$tdoubling$ci$cex,
        font = control$tdoubling$ci$font
      )
      ey <- add_lines_to_user(0.15, add_height_to_user(cih, ciy, log), log)
      eh <- strheight("",
        units = "user",
        cex = control$tdoubling$estimate$cex,
        font = control$tdoubling$estimate$font
      )
      e <- log(2) / exp(cache_log_r$estimate)
      l <- log(2) / exp(cache_log_r$upper)
      u <- log(2) / exp(cache_log_r$lower)
      n <- length(wl)
      get_el <- function(el) sapply(control$tdoubling[c("estimate", "ci")], `[[`, el)
      text(
        x = rep((t1 + t2) / 2, times = 2L),
        y = rep(c(ey, ciy), each = n),
        labels = c(sprintf("%.1f", e), sprintf("(%.1f, %.1f)", l, u)),
        cex = rep(cex * get_el("cex"), each = n),
        col = rep(get_el("col"), each = n),
        font = rep(get_el("font"), each = n),
        adj = c(0.5, 0),
        xpd = TRUE
      )
      ciy_cap <- add_lines_to_user(0.5, add_height_to_user(eh, ey, log), log)
      ey_cap <- add_lines_to_user(0.15, add_height_to_user(cih, ciy_cap, log), log)
      ew <- strwidth("estimate",
        units = "user",
        cex = control$tdoubling$estimate$cex,
        font = control$tdoubling$estimate$font
      )
      text(
        x = rep_len(usr[2L] - 0.5 * ew, 2L),
        y = c(ey_cap, ciy_cap),
        labels = c("estimate", "(95% CI)"),
        cex = cex * get_el("cex"),
        col = get_el("col"),
        font = get_el("font"),
        adj = c(0.5, 0),
        xpd = TRUE
      )
      text(
        x = usr[2L],
        y = add_lines_to_user(0.25, add_height_to_user(eh, ey_cap, log), log),
        labels = "initial doubling time, days:",
        cex = cex * control$tdoubling$caption$cex,
        col = control$tdoubling$caption$col,
        font = control$tdoubling$caption$font,
        adj = c(1, 0),
        xpd = TRUE
      )
    }

    ## Box
    if (!is.null(control$box)) {
      do.call(box, control$box)
    }

    ## Axis (x)
    if (!is.null(control$axis$x)) {
      if (time_as == "Date") {
        daxis(
          left = usr[1L],
          right = usr[2L],
          origin = origin + shift + 1,
          minor = control$axis$x$minor,
          major = control$axis$x$major
        )
      } else {
        l <- list(side = 1)
        do.call(axis, c(l, control$axis$x))
      }
    }

    ## Axis (y)
    if (!is.null(control$axis$y)) {
      l1 <- list(side = 2, las = 1)
      l2 <- control$axis$y
      l1$at <- axTicks(side = 2)
      if (type != "rt1" && (max(l1$at) >= 1e05)) {
        l1$labels <- get_yax_labels(l1$at)
        l2$cex.axis <- min(l2$cex.axis, get_yax_cex(l1$labels, mex = 3.5 - l2$mgp[2L]))
      }
      do.call(axis, c(l1, l2))
      if (type == "rt1" && usr[4L] > 0) {
        axis(
          side = 2,
          line = 4.5,
          at = c(0, usr[4L]),
          labels = c("", ""),
          lwd = l2$lwd,
          col = l2$col,
          lwd.ticks = 0,
          xpd = TRUE
        )
        tdoubling <- c(1:5, 10, 50, 100)
        l1$at <- log(2) / tdoubling
        l1$labels <- tdoubling
        ## FIXME: why `mgp[4] = 4 + foo`, not `mgp[2] = 4.5 + foo`?
        l2$mgp <- c(3, 4 + l2$mgp[2L], 4.5)
        l2$lwd <- 0
        do.call(axis, c(l1, l2))
      }

      ## Axis title (x)
      if (!is.null(control$title$xlab)) {
        l <- list(xlab = xlab, line = 2.5)
        do.call(title, c(l, control$title$xlab))
      }

      ## Axis title (y)
      if (!is.null(control$title$ylab)) {
        if (type == "rt1") {
          mlw <- max(strwidth(axTicks(side = 2),
            units = "inches",
            cex = control$axis$y$cex.axis,
            font = control$axis$y$font.axis
          ))
          line0 <- 0.5 + mlw / csi + control$axis$y$mgp[2L]
        } else {
          line0 <- 4
        }
        l <- list(ylab = ylab, line = line0)
        do.call(title, c(l, control$title$ylab))
        if (type == "rt1" && usr[4L] > 0) {
          mlw <- max(strwidth(tdoubling,
            units = "inches",
            cex = control$axis$y$cex.axis,
            font = control$axis$y$font.axis
          ))
          tw <- strwidth(ylab_outer,
            units = "inches",
            cex = control$title$ylab$cex.lab,
            font = control$title$ylab$font.lab
          )
          ## `tw` relative to 0.8 times the distance from 0 to `usr[4]`
          rtw <- (tw * diff(usr[3:4]) / pin[2L]) / (0.8 * (usr[4L] - max(0, usr[3L])))
          text(
            ## FIXME: why 7.5, not 5 = 4.5 + 0.5?
            x = usr[1L] - cxy[1L] * (7.5 + mlw / csi + control$axis$y$mgp[2L]),
            y = mean(c(max(0, usr[3L]), usr[4L])),
            labels = ylab_outer,
            adj = c(0.5, 0),
            srt = 90,
            xpd = NA,
            cex = cex * control$title$ylab$cex.lab / max(1, rtw),
            col = control$title$ylab$col.lab,
            font = control$title$ylab$font.lab
          )
        }
      }
    }

    ## Plot (sub)title
    if (!is.null(control$title$main)) {
      if (show_tdoubling) {
        line0 <- user_to_lines(add_lines_to_user(0.5, add_height_to_user(eh, ey, log), log), log)
      } else {
        line0 <- 0.5
      }
      if (!is.null(sub) && !is.null(control$title$sub)) {
        names(control$title$sub) <- base::sub("\\.sub$", ".main", names(control$title$sub))
        th <- strheight(sub[k],
          units = "inches",
          cex = control$title$sub$cex.main,
          font = control$title$sub$font.main
        )
        l <- list(main = sub[k], line = line0)
        do.call(title, c(l, control$title$sub))
        line0 <- line0 + th / csi + 0.25
      }
      l <- list(main = main[k], line = line0)
      do.call(title, c(l, control$title$main))
    }

    ## Legend
    if (show_legend) {
      lx <- usr[2L] + 0.02 * diff(usr[1:2])
      ly <- inv_log(usr[4L] - 0.02 * diff(usr[3:4]), log)
      cond <- all(data$dt[data$pty == pty[1L]] == m, na.rm = TRUE)
      lexpr <- parse(
        text = sprintf("'%s,' ~ Delta * 't %s %g day%s'",
          c("obs", "obs", "obs", "pred"),
          c((if (cond) "=" else "~"), "<", ">", "="),
          m,
          if (m > 1) "s" else ""
        )
      )
      lind <- c(pty %in% levels(data$pty), TRUE)
      get_el <- function(el) sapply(control$special[pty], `[[`, el)
      legend(
        x = lx,
        y = ly,
        legend = lexpr[lind],
        xpd = NA,
        bty = "n",
        cex = 0.7,
        seg.len = 1,
        pch = c(get_el("pch"), NA)[lind],
        pt.bg = c(get_el("bg"), NA)[lind],
        lty = c(NA, NA, NA, control$lines$lty)[lind],
        lwd = c(NA, NA, NA, control$lines$lwd)[lind],
        col = c(get_el("col"), control$lines$col)[lind]
      )
    }
  }

  invisible(NULL)
}

#' @import graphics
#' @importFrom grDevices colorRamp rgb
do_heat_plot <- function(cache, origin, time_as, per_plot, control,
                         xlim, main, xlab, ylab, ylab_outer, plab) {
  ## Axis limits
  tr <- range(cache$time, na.rm = FALSE)
  shift <- tr[1L]
  if (is.null(xlim)) {
    xlim <- tr - shift
  } else {
    if (inherits(xlim, "Date")) {
      xlim <- julian(xlim, origin = origin + shift)
    }
  }
  ylim <- c(0, 1)

  ## Titles
  if (is.null(main)) {
    main <- "Per capita growth rate, by time series"
  }
  if (is.null(xlab)) {
    xlab <- switch(time_as,
      Date = "",
      sprintf("number of days since %s", origin + shift)
    )
  }
  if (is.null(ylab)) {
    ylab <- "growth\nrate,\nper day"
  }
  if (is.null(ylab_outer)) {
    ylab_outer <- "doubling\ntime,\ndays"
  }
  if (is.null(plab)) {
    plab <- levels(cache$ts)
  }

  ## Heat map color palette
  pal <- do.call(colorRamp, control$colorRamp)
  ips <- control$ips

  op <- par(oma = c(2.5, 0, 3.5, 0))
  on.exit(par(op))

  ### Loop over plots =====================================

  cache$time <- cache$time - shift
  cache$estimate <- exp(cache$estimate)
  cache_split <- split(cache, cache$ts)
  tsl <- levels(cache$ts)
  N <- length(tsl)
  M <- max(cache$estimate, na.rm = TRUE)

  K <- 0L
  while (K < N) {
    L <- c(seq_len(per_plot), rep_len(per_plot + 1L, per_plot))
    dim(L) <- c(per_plot, 2L)
    layout(L, widths = c(par("din")[1L] - 1.2, 1.2))
    par(mar = c(0.5 * ips, 1.5, 0.5 * ips, 0.5))

    ### Loop over panels ==================================

    for (k in K + seq_len(min(per_plot, N - K))) {
      plot.new()
      plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i")

      cxy <- par("cxy")
      pin <- par("pin")
      usr <- par("usr")

      d <- cache_split[[k]]

      ## Panel background
      polygon(
        x = xlim[c(1L, 2L, 2L, 1L)],
        y = ylim[c(1L, 1L, 2L, 2L)],
        border = NA,
        col = 1
      )

      ## Panel pixels
      for (i in seq_len(nrow(d))) {
        polygon(
          x = (d$time[i] + c(-0.5, 0.5))[c(1L, 2L, 2L, 1L)],
          y = ylim[c(1L, 1L, 2L, 2L)],
          border = NA,
          col = rgb(pal(d$estimate[i] / M), maxColorValue = 255)
        )
      }

      ## Panel title background
      px <- diff(usr[1:2]) * (0.075 * pin[2L] / pin[1L])
      py <- 0.075 * diff(usr[3:4])
      tw <- strwidth(plab[k], cex = 1, font = 2)
      th <- strheight(plab[k], cex = 1, font = 2)
      rect(
        xleft = usr[1L],
        ybottom = usr[4L] - th - 2 * px,
        xright = usr[1L] + tw + 2 * py,
        ytop = usr[4L],
        border = NA,
        col = "#00000080"
      )

      ## Panel title
      text(
        x = usr[1L] + px,
        y = usr[4L] - py,
        labels = plab[k],
        adj = c(0, 1)
      )

      ## Plot title
      if (k == K + 1L) {
        title(
          main = main,
          line = 0.5,
          adj = 0,
          cex.main = 1.3,
          font = 2,
          xpd = NA
        )
      }
    }

    ## Axis (x)
    axis(
      side = 1,
      at = usr[1:2],
      labels = c("", ""),
      lwd.ticks = 0
    )
    if (time_as == "Date") {
      daxis(
        left = xlim[1L],
        right = xlim[2L],
        origin = origin + shift + 1,
        minor = list(mgp = c(3, 0.25, 0), tcl = -0.2, cex.axis = 1,   gap.axis = 0, xpd = TRUE),
        major = list(mgp = c(3, 1.25, 0), tcl = 0,    cex.axis = 1.3, gap.axis = 0, xpd = TRUE)
      )
    } else {
      axis(side = 1, mgp = c(3, 0.7, 0), gap.axis = 0, xpd = TRUE)
    }

    ## Axis title (x)
    title(xlab = xlab, line = 2.5)

    ## Skip empty panels to get to the last
    for (i in seq_len((per_plot - (K %% per_plot)) %% per_plot)) {
      plot.new()
    }

    ## Color scale
    par(mar = c(0.5 * ips, 1, 0.5 * ips, 7))
    plot.new()
    dy2 <- 0.0025
    plot.window(
      xlim = c(0, 1),
      ylim = c(-dy2, 1 + dy2) * M,
      xaxs = "i",
      yaxs = "i"
    )
    for (y in seq.int(0, 1, by = 2 * dy2)) {
      polygon(
        x = c(0, 1)[c(1L, 2L, 2L, 1L)],
        y = ((y + c(-dy2, dy2)) * M)[c(1L, 1L, 2L, 2L)],
        border = NA,
        col = rgb(pal(y), maxColorValue = 255)
      )
    }

    ## Axis (y, inner)
    axis(
      side = 4,
      at = usr[3:4],
      labels = c("", ""),
      lwd.ticks = 0
    )
    axis(
      side = 4,
      mgp = c(3, 0.7, 0),
      las = 1,
      lwd = 0,
      lwd.ticks = 1
    )

    ## Axis (y, outer)
    tdoubling <- c(1:5, 10, 50, 100)
    axis(
      side = 4,
      at = usr[3:4],
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

    ## Axis title (y, inner)
    text(
      x = usr[2L] + 0 * cxy[1L],
      y = usr[4L] + 0.5 * cxy[2L],
      labels = ylab,
      cex = 1,
      adj = c(0, 0),
      xpd = NA
    )

    ## Axis title (y, outer)
    text(
      x = usr[2L] + 3.5 * cxy[1L],
      y = usr[4L] + 0.5 * cxy[2L],
      labels = ylab_outer,
      cex = 1,
      adj = c(0, 0),
      xpd = NA
    )

    K <- K + per_plot
  }

  invisible(NULL)
}
