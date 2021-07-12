#' Plot nonlinear mixed effects models of epidemic growth
#'
#' Methods for plotting \code{"\link{egf}"} objects.
#'
#' @param x
#'   An \code{"\link{egf}"} object.
#' @param type
#'   A \link{character} string indicating a type of plot. The options are:
#'   interval incidence (\code{"interval"}),
#'   cumulative incidence (\code{"cumulative"}),
#'   per capita growth rate as curve (\code{"rt1"}), and
#'   per capita growth rate as heat map (\code{"rt2"}).
#'   \code{"rt2"} displays \code{panels_per_plot} time series per plot.
#'   The rest display one time series per plot.
#' @param subset
#'   An expression to be evaluated in the combined model frame
#'   (see \code{\link{make_combined}}). Must evaluate to
#'   a \link{logical} vector indexing rows of the data frame,
#'   and thus fitting windows. Only indexed fitting windows are plotted.
#'   The default (\code{\link{NULL}}) is to plot all fitting windows.
#' @param order
#'   An expression to be evaluated in the combined model frame
#'   (see \code{\link{make_combined}}), typically a call to
#'   \code{\link{order}}, determining the order in which time series
#'   are plotted. The default (\code{\link{NULL}}) is equivalent to
#'   \code{\link{seq_len}(\link{nrow}(combined))}.
#' @param cache
#'   An \code{"egf_plot_cache"} object returned by a previous evaluation
#'   of \code{plot.egf(x)}. Fitted and predicted values and standard errors
#'   stored in \code{cache} are reused if possible to avoid recomputation.
#' @param do_plot
#'   A \link{logical} flag. If \code{FALSE}, then nothing is plotted.
#'   Useful when only the returned \code{"egf_plot_cache"} object is desired.
#' @param time_as
#'   A \link{character} string indicating how time is displayed on the
#'   bottom axis. The options are: as a calendar (\code{"Date"}) and
#'   as a number of days since the earliest time point (\code{"numeric"}).
#' @param dt
#'   A positive number indicating an observation interval in days.
#'   Predicted curves are evaluated on grids with this spacing.
#'   When \code{type = "interval"}, counts observed over a shorter
#'   or longer interval \code{dt0} are scaled by a factor of
#'   \code{dt / dt0} so that their scale matches that of the curves.
#'   These points can be highlighted via \code{control}.
#' @param log
#'   A \link{logical} flag. If \code{TRUE}, then the dependent variable
#'   is plotted on a logarithmic scale.
#'   [\code{type != "rt1"} only.]
#' @param show_predict
#'   An integer flag: 2 is to draw predicted curves with confidence bands,
#'   1 is draw predicted curves only, 0 is to draw neither.
#'   \link[=logical]{Logical} values are coerced to integer.
#'   [\code{type != "rt2"} only.]
#' @param show_tdoubling
#'   An integer flag: 2 is to print initial doubling time estimates in
#'   the top margin with confidence intervals, 1 is to print estimates
#'   only, 0 is to print neither. \link[=logical]{Logical} values are
#'   coerced to integer. Supported only if \code{x$model$curve} is
#'   \code{"exponential"}, \code{"logistic"}, or \code{"richards"}.
#'   [\code{type != "rt2"} only.]
#' @param show_legend
#'   A \link{logical} flag. If \code{TRUE}, then a legend is displayed
#'   in the right margin.
#'   [\code{type = "interval"} only.]
#' @param level
#'   A number in the interval (0,1). This is the confidence level used
#'   when \code{show_predict = 2} or \code{show_tdoubling = 2}.
#'   [\code{type != "rt2"} only.]
#' @param panels_per_plot
#'   A positive integer giving the number of panels (time series)
#'   displayed in one plot.
#'   [\code{type = "rt2"} only.]
#' @param inter_panel_space
#'   A non-negative number giving the space between panels
#'   in multipanel plots as a number of margin lines.
#'   [\code{type = "rt2"} only.]
#' @param control
#'   An \code{"\link{egf_plot_control}"} object controlling the appearance
#'   of almost all plot elements.
#' @param xlim,ylim
#'   \link[=numeric]{Numeric} vectors of length 2 specifying axis limits,
#'   which are recycled for all plots. If \code{time_as = "Date"}, then
#'   \code{xlim} can instead be a \link{Date} vector or a \link{character}
#'   vector coercible to Date via \code{\link{as.Date}(xlim)}. \code{ylim}
#'   is unused by \code{type = "rt2"}.
#' @param main,sub,xlab,ylab,ylab_outer,plab
#'   \link[=character]{Character} strings or expressions used
#'   to generate plot (\code{main}, \code{sub}), axis (\code{xlab},
#'   \code{ylab}, \code{ylab_outer}), and panel (\code{plab}) titles.
#'   \code{main}, \code{xlab}, and \code{ylab} are supported for all
#'   values of \code{type}.
#'   \code{sub} is unused by \code{type = "rt2"}.
#'   \code{plab} is used by \code{type = "rt2"} only.
#'   \code{ylab_outer} is used by \code{type = "rt[12]"} only.
#'   When \code{type != "rt2"}, \code{main} and \code{sub} are evaluated
#'   in the combined model frame (see \code{\link{make_combined}})
#'   in order to generate unique (sub)titles for each plot.
#'   When \code{type = "rt2"}, \code{plab} is evaluated similarly
#'   in order to generate unique titles for each panel.
#'   \code{\link{plotmath}} expressions are not supported
#'   for \code{main}, \code{sub}, and \code{plab} in these cases.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A \link[=data.frame]{data.frame} inheriting from \link{class}
#' \code{"egf_plot_cache"}. If argument \code{cache} was supplied
#' in the function call, then this data frame is the result of
#' augmenting \code{cache} with new computations.
#'
#' @details
#' Computation of fitted and predicted values and standard errors
#' is performed before any plots are created. To avoid waste of
#' computation time, cached computations are returned \emph{even if}
#' an error is thrown during plotting. The cache will be temporarily
#' available via \code{\link{.Last.value}}. To ensure that the cache
#' is permanently available, assign the result of the call to
#' \code{\link{plot}} to a name: \code{cache <- plot(x, \dots)}.
#'
#' Caching functionality must be used with care, as mismatch
#' between \code{x} and \code{cache} will \emph{not} be detected.
#' Constructions such as \code{plot(y, cache = plot(x, \dots), \dots)}
#' should \emph{not} be expected to produce correct results.
#'
#' See topic \code{\link{nse}} for details on nonstandard evaluation
#' of \code{subset} and \code{order}.
#'
#' @export
#' @importFrom stats fitted predict complete.cases
plot.egf <- function(x,
                     type = c("interval", "cumulative", "rt1", "rt2"),
                     subset = NULL,
                     order = NULL,
                     cache = NULL,
                     do_plot = TRUE,
                     time_as = c("Date", "numeric"),
                     dt = 1,
                     log = TRUE,
                     show_predict = TRUE,
                     show_tdoubling = TRUE,
                     show_legend = FALSE,
                     level = 0.95,
                     panels_per_plot = min(6L, nlevels(x$frame$ts)),
                     inter_panel_space = 0.25,
                     control = egf_plot_control(),
                     xlim = NULL,
                     ylim = NULL,
                     main = NULL,
                     sub = NULL,
                     xlab = NULL,
                     ylab = NULL,
                     ylab_outer = NULL,
                     plab = NULL,
                     ...) {
  ## FIXME: `order` _actually_ orders `levels(window)`,
  ## not `levels(ts)`, leading to some indexing gore ...

  type <- match.arg(type)
  stop_if_not_true_false(do_plot)
  stop_if_not_number(dt, "positive")
  if (type == "rt2") {
    show_predict <- 1L
    show_tdoubling <- 0L
  } else {
    stop_if_not_true_false(show_predict, allow_numeric = TRUE)
    show_predict <- min(2L, max(0L, as.integer(show_predict))) # coercion to `0:2`
    if (x$model$curve %in% c("exponential", "logistic", "richards")) {
      stop_if_not_true_false(show_tdoubling, allow_numeric = TRUE)
      show_tdoubling <- min(2L, max(0L, as.integer(show_tdoubling))) # coercion to `0:2`
    } else {
      show_tdoubling <- 0L
    }
    if (any(c(show_predict, show_tdoubling) == 2L)) {
      stop_if_not_number_in_interval(level, 0, 1, "()")
    }
  }

  combined <- make_combined(x)
  subset <- subset_to_index(substitute(subset), data = combined, enclos = parent.frame())
  stop_if_not(
    sum(subset) > 0L,
    m = "`subset` must index at least one fitting window."
  )

  ## This code _could_ be run _after_ augmenting `cache`,
  ## but running it _before_ avoids waste of computation
  ## time in the event of errors. This really only matters
  ## for users who have not assigned `plot(x)`.
  if (do_plot) {
    time_as <- match.arg(time_as)

    if (type == "interval") {
      stop_if_not_true_false(show_legend)
    } else {
      show_legend <- FALSE
    }
    if (type == "rt1") {
      log <- FALSE
    } else {
      stop_if_not_true_false(log)
    }
    if (type == "rt2") {
      stop_if_not_integer(panels_per_plot, "positive")
      stop_if_not_number(inter_panel_space, "nonnegative")
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

    stop_if_not(
      inherits(control, "egf_plot_control"),
      m = "`control` must inherit from class \"egf_plot_control\". See `?egf_plot_control`."
    )

    order <- order_to_index(substitute(order), data = combined, enclos = parent.frame())
    if (type == "rt2") {
      plab <- label_to_character(substitute(plab), data = combined, enclos = parent.frame())
    } else {
      main <- label_to_character(substitute(main), data = combined, enclos = parent.frame())
      sub <- label_to_character(substitute(sub), data = combined, enclos = parent.frame())
    }

    frame <- x$frame
    wl <- as.character(x$endpoints$window) # same as levels
    tsl <- as.character(x$endpoints$ts) # same as levels but with duplicates

    wl_subset <- wl[subset]
    frame$window <- factor(frame$window, levels = wl_subset)

    tsl_subset <- intersect(tsl[order], frame$ts[!is.na(frame$window)])
    frame$ts <- factor(frame$ts, levels = tsl_subset)

    if (type == "cumulative") {
      keep <- !c(tapply(frame$x, frame$ts, function(x) anyNA(x[-1L])))
      if (any(!keep)) {
        warning(
          wrap(
            "Missing values preventing calculation of cumulative incidence. ",
            "These time series will not be displayed:"
          ),
          "\n\n", paste0("  ", tsl_subset[!keep], collapse = "\n")
        )
        tsl_subset <- tsl_subset[keep]
        frame$ts <- factor(frame$ts, levels = tsl_subset)
        subset <- subset & (tsl %in% tsl_subset)
        wl_subset <- wl[subset]
        frame$window <- factor(frame$window, levels = wl_subset)
      }
      if (all(!keep)) {
        stop("Nothing left to do ...")
      }
    }
  }

  ## If necessary, initialize cache
  if (is.null(cache)) {
    cache <- data.frame(
      var = character(0L),
      ts = character(0L),
      window = character(0L),
      time = numeric(0L),
      estimate = numeric(0L),
      se = numeric(0L),
      stringsAsFactors = TRUE
    )
  } else {
    stop_if_not(
      inherits(cache, "egf_plot_cache"),
      m = "`cache` must inherit from class \"egf_plot_cache\"."
    )
  }
  nc <- names(cache)

  ## If necessary, augment `cache` with predicted values
  ## of whatever is being plotted and standard errors
  if (show_predict > 0L) {
    what <- switch(type, interval = "log_int_inc", cumulative = "log_cum_inc", "log_rt")
    ok <- cache$var == what & !(show_predict == 2L & is.na(cache$se))
    need <- setdiff(x$endpoints$window[subset], cache$window[ok])
    if (length(need) > 0L) {
      w <- match(need, x$endpoints$window, 0L)
      time_split <- Map(seq.int,
        from = x$endpoints$start[w],
        to = x$endpoints$end[w],
        by = dt
      )
      pd <- predict(x,
        what = what,
        time = unlist(time_split, FALSE, FALSE),
        window = rep.int(x$endpoints$window[w], lengths(time_split)),
        log = TRUE,
        se = (show_predict == 2L)
      )
      if (show_predict == 1L) {
        pd$se <- NA_real_
      }
      cache <- rbind(cache, pd[nc])
    }
  }

  ## If necessary, augment `cache` with fitted values of `log(r)`
  ## and standard errors
  if (show_tdoubling > 0L) {
    ok <- cache$var == "log(r)" & !(show_tdoubling == 2L & is.na(cache$se))
    need <- setdiff(x$endpoints$window[subset], cache$window[ok])
    if (length(need) > 0L) {
      ft <- fitted(x,
        par = "log(r)",
        link = TRUE,
        se = (show_tdoubling == 2L),
        .subset = (x$endpoints$window %in% need)
      )
      if (show_tdoubling == 1L) {
        ft$se <- NA_real_
      }
      ft$time <- NA_real_
      names(ft)[match("par", names(ft), 0L)] <- "var"
      cache <- rbind(cache, ft[nc])
    }
  }

  ## Clean up
  o <- do.call(base::order, unname(cache[nc[c(1:4, 6L)]]))
  cache <- cache[o, , drop = FALSE]
  i <- !duplicated(cache[nc[1:4]])
  cache <- cache[i, , drop = FALSE]
  row.names(cache) <- NULL
  class(cache) <- c("egf_plot_cache", "data.frame")

  ## If not plotting, then return
  if (!do_plot) {
    return(invisible(cache))
  }

  ## If plotting, then create an instruction to return
  ## `cache` if the low level plot function throws an error
  cache_bak <- cache
  on.exit({
    message("Augmented `cache` returned despite error ...")
    return(invisible(cache_bak))
  })

  ## Extract only those rows of `cache` needed by the low level plot function
  cache[nc[1:3]] <- Map(factor,
    x = cache[nc[1:3]],
    levels = list(c(what, if (show_tdoubling > 0L) "log(r)"), tsl_subset, wl_subset)
  )
  i <- complete.cases(cache[nc[1:3]])
  cache <- cache[i, , drop = FALSE]

  ## Compute confidence intervals
  i <- (show_predict == 2L & cache$var == what) | (show_tdoubling == 2L & cache$var == "log(r)")
  cache[c("lower", "upper")] <- NA_real_
  cache[i, c("lower", "upper")] <- do_wald(
    estimate = cache$estimate[i],
    se = cache$se[i],
    level = level
  )

  ## Subset and order plot/axis/panel titles
  ## so that they correspond elementwise to `tsl_subset`
  subset <- order[order %in% subset]
  m <- match(tsl_subset, tsl[subset], 0L)

  origin <- attr(x$frame, "origin")

  if (type == "rt2") {
    do_heat_plot(
      cache = cache,
      time_as = time_as,
      dt = dt,
      origin = origin,
      log = log,
      panels_per_plot = panels_per_plot,
      inter_panel_space = inter_panel_space,
      control = control,
      xlim = xlim,
      main = main,
      sub = sub,
      xlab = xlab,
      ylab = ylab,
      ylab_outer = ylab_outer,
      plab = plab[subset][m]
    )
  } else {
    do_curve_plot(
      frame = frame,
      cache = cache,
      type = type,
      time_as = time_as,
      dt = dt,
      origin = origin,
      log = log,
      curve = x$model$curve,
      show_predict = show_predict,
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

  ## Discard exit instructions if low level plot function runs without stopping
  on.exit()
  invisible(cache)
}

#' @import graphics
do_curve_plot <- function(frame, cache, type, time_as,
                          dt, origin, log, curve,
                          show_predict, show_tdoubling, show_legend,
                          level, control, xlim, ylim,
                          main, sub, xlab, ylab, ylab_outer) {
  ### Set up ===================================================================

  frame_split <- split(frame, frame$ts)
  cache_split <- split(cache, cache$ts)
  N <- nlevels(frame$ts)
  formula <- as.formula(call("~", as.name(type), quote(time)))
  var <- switch(type, interval = "log_int_inc", cumulative = "log_cum_inc", "log_rt")

  xlim_bak <- xlim
  ylim_bak <- ylim
  xlab_bak <- xlab

  ## Plot titles
  if (is.null(main)) {
    if (curve %in% c("gompertz", "richards")) {
      substr(curve, 1L, 1L) <- toupper(substr(curve, 1L, 1L))
    }
    main <- rep_len(sprintf("Fitted %s model", curve), N)
  }

  ## Plot subtitles
  if (is.null(sub)) {
    sub <- names(frame_split)
  }

  ## Axis title (y)
  if (is.null(ylab)) {
    if (type == "rt1") {
      ylab <- "growth rate, per day"
    } else {
      ylab <- paste(type, "incidence")
    }
  }
  if (is.null(ylab_outer)) {
    ylab_outer <- "doubling time, days"
  }

  ## Utilities
  win_to_lines <- function(win) {
    diff(grconvertX(c(0, win), "inches", "lines"))
  }
  hin_to_lines <- function(hin) {
    diff(grconvertY(c(0, hin), "inches", "lines"))
  }
  yu_to_y <- if (log) function(yu) 10^yu else identity

  ## Hack for `points` loop
  control$points_basic <- control$points

  ## Graphical parameters
  if (type == "rt1") {
    ## Add space for secondary axis in left margin if necessary
    mar <- c(3.5, 4 + (is.null(ylim) || ylim[2L] > 0) * 4, 4, 1) + 0.1
  } else {
    ## Add space for legend in right margin if necessary
    mar <- c(3.5, 5, 4, 1 + 5.5 * show_legend) + 0.1
  }
  op <- par(mar = mar, xaxs = "i", yaxs = "i")
  on.exit(par(op))


  ### Loop over plots ==========================================================

  for (k in seq_len(N)) {
    ### Set up for plot --------------------------------------------------------

    frame <- droplevels(frame_split[[k]])
    shift <- min(frame$time)
    n <- nlevels(frame$window)

    data <- data.frame(
      time = frame$time - shift,
      dt = c(NA, diff(frame$time))
    )
    data[[type]] <- switch(type,
      interval = c(NA, frame$x[-1L]) * dt / data$dt,
      cumulative = c(0, cumsum(frame$x[-1L])),
      rt1 = c(NA, diff(base::log(frame$x))) / data$dt
    )
    data[[type]][!is.finite(data[[type]])] <- NA

    t1 <- c(tapply(data$time, frame$window, min))
    t2 <- c(tapply(data$time, frame$window, max))

    cache <- droplevels(cache_split[[k]])
    cache_log_r <- cache[cache$var == "log(r)", , drop = FALSE]
    cache_predict <- cache[cache$var == var, , drop = FALSE]
    cache_predict$time <- cache_predict$time - shift
    cache_predict_split <- split(cache_predict, cache_predict$window)

    ## Axis limits (x)
    if (is.null(xlim_bak)) {
      xlim <- c(0, max(data$time) * 1.01)
    } else if (inherits(xlim_bak, "Date")) {
      xlim <- julian(xlim_bak, origin = origin + shift)
    } else {
      xlim <- xlim_bak
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
    } else {
      ylim <- ylim_bak
    }

    ## Axis title (x)
    if (is.null(xlab_bak)) {
      xlab <- switch(time_as,
        Date = "",
        sprintf("number of days since %s", origin + shift)
      )
    } else {
      xlab <- xlab_bak
    }

    ## Point styles
    if (type == "interval") {
      data$pty <- factor(sign(data$dt - dt), levels = c(0, -1, 1), labels = c("basic", "short", "long"))
    } else {
      data$pty <- factor("basic")
    }


    ### Plot -------------------------------------------------------------------

    plot.new()
    plot.window(xlim = xlim, ylim = ylim, log = if (log) "y" else "")
    usr <- par("usr")
    pin <- par("pin")
    cex <- par("cex")

    ## Fitting windows
    if (!is.null(ct <- control$rect)) {
      args <- list(
        xleft   = t1,
        ybottom = yu_to_y(usr[3L]),
        xright  = t2,
        ytop    = yu_to_y(usr[4L])
      )
      do.call(rect, c(args, ct))
    }

    ## Observed data
    for (s in levels(data$pty)) {
      if (!is.null(ct <- control[[sprintf("points_%s", s)]])) {
        args <- list(formula = formula, data = data, subset = (data$pty == s))
        do.call(points, c(args, ct))
      }
    }

    ## Confidence bands on predicted curves
    if (show_predict == 2L && !is.null(ct <- control$polygon)) {
      for (pd in cache_predict_split) {
        pd$lower <- exp(pd$lower)
        pd$upper <- exp(pd$upper)
        if (type == "cumulative") {
          c0 <- data[[type]][match(pd$time[1L], data$time, 0L)]
          pd$lower <- c0 + pd$lower
          pd$upper <- c0 + pd$upper
        }
        args <- list(
          x = c(pd$time, rev(pd$time)),
          y = c(pd$lower, rev(pd$upper))
        )
        do.call(polygon, c(args, ct))
      }
    }

    ## Predicted curves
    if (show_predict > 0L && !is.null(ct <- control$lines)) {
      for (pd in cache_predict_split) {
        pd$estimate <- exp(pd$estimate)
        if (type == "cumulative") {
          c0 <- data[[type]][match(pd$time[1L], data$time, 0L)]
          pd$estimate <- c0 + pd$estimate
        }
        args <- list(formula = estimate ~ time, data = pd)
        do.call(lines, c(args, ct))
      }
    }

    ## Asymptotes
    if (type == "rt1") {
      if (!is.null(ct <- control$abline)) {
        args <- list(h = 0)
        do.call(abline, c(args, ct))
      }
      if (!is.null(ct <- control$segments)) {
        r <- exp(cache_log_r$estimate)
        args <- list(x0 = t1, y0 = r, x1 = t2, y1 = r)
        do.call(segments, c(args, ct))
      }
    }

    ## Box
    if (!is.null(ct <- control$box)) {
      do.call(box, ct)
    }

    if (time_as == "Date") {
      ## Axis (x)
      Daxis(
        origin = origin + shift + 1,
        minor = control$axis_x_Date_minor,
        major = control$axis_x_Date_major,
        show_minor = !is.null(control$axis_x_Date_minor),
        show_major = !is.null(control$axis_x_Date_major)
      )
    } else {
      ## Axis (x)
      if (!is.null(ct <- control$axis_x_numeric)) {
        args <- list(side = 1)
        do.call(baxis, c(args, ct1))
      }

      ## Axis title (x)
      if (!is.null(ct <- control$title_xlab)) {
        args <- list(xlab = xlab, line = 2.5)
        do.call(title, c(args, ct))
      }
    }

    if (!is.null(ct1 <- control$axis_y)) {
      ## Axis (y)
      args <- list(side = 2, at = axTicks(side = 2), las = 1)
      if (type != "rt1" && max(args$at) >= 1e05) {
        args$labels <- get_yax_labels(args$at)
        yac <- get_yax_cex(args$labels,
          mex = 3.5 - ct1$mgp[2L],
          font = ct1$font.axis,
          family = ct1$family
        )
        ct1$cex.axis <- min(yac, ct1$cex.axis)
      }
      do.call(baxis, c(args, ct1))

      if (type == "rt1" && usr[4L] > 0) {
        tdoubling <- c(1:5, 10, 20, 50, 100)
        args <- list(
          side = 2,
          a = 0,
          b = usr[4L],
          at = log(2) / tdoubling,
          labels = tdoubling,
          las = 1
        )
        ct1_outer <- ct1
        ct1_outer$mgp <- ct1_outer$mgp + 4
        do.call(baxis, c(args, ct1_outer))
      }

      ## Axis title (y)
      if (!is.null(ct2 <- control$title_ylab)) {
        if (type == "rt1") {
          win_tick_labels <- strwidth(axTicks(side = 2), units = "inches", cex = ct1$cex.axis, font = ct1$font.axis, family = ct1$family)
          line <- 0.5 + win_to_lines(max(win_tick_labels)) + ct1$mgp[2L]
        } else {
          line <- 4
        }
        args <- list(ylab = ylab, line = line)
        do.call(title, c(args, ct2))

        if (type == "rt1" && usr[4L] > 0) {
          win_tick_labels <- strwidth(tdoubling, units = "inches", cex = ct1$cex.axis, font = ct1$font.axis, family = ct1$family)
          win_axis_title <- strwidth(ylab_outer, units = "inches", cex = ct2$cex.lab, font = ct2$font.lab, family = ct2$family)
          line_ylab_outer <- 0.5 + win_to_lines(max(win_tick_labels)) + ct1_outer$mgp[2L]
          rho <- (usr[4L] - max(0, usr[3L])) /
            (win_axis_title * (usr[4L] - usr[3L]) / pin[2L])
          mtext(ylab_outer,
            side = 2,
            line = line_ylab_outer,
            at = (max(0, usr[3L]) + usr[4L]) / 2,
            xpd = NA,
            col = ct2$col.lab,
            cex = ct2$cex.lab * min(1, 0.8 * rho),
            font = ct2$font.lab,
            family = ct2$family
          )
        }
      }
    }

    ## Initial doubling times
    if (show_tdoubling > 0L) {
      s <- c("estimate", "upper", "lower")
      elu <- base::log(2) / exp(cache_log_r[s])
      names(elu) <- s[c(1L, 3L, 2L)]

      s <- c("caption", "estimate", "ci")
      ct <- control[sprintf("mtext_tdoubling_%s", s)]
      names(ct) <- s

      show_caption <- !is.null(ct$caption)

      hin_ci <- strheight("", units = "inches", cex = ct$ci$cex, font = ct$ci$font, family = ct$ci$family)
      hin_e <- strheight("", units = "inches", cex = ct$estimate$cex, font = ct$estimate$font, family = ct$estimate$family)

      line_ci     <- 0.25
      line_e      <- line_ci    + hin_to_lines(hin_ci) + 0.15
      line_lg_ci  <- line_e     + hin_to_lines(hin_e)  + 0.5
      line_lg_e   <- line_lg_ci + hin_to_lines(hin_ci) + 0.15
      line_lg_cap <- line_lg_e  + hin_to_lines(hin_e)  + 0.25

      adj <- ct$caption$adj
      if (is.null(adj)) {
        adj <- 1
      }

      wu_lg_cap <- strwidth("initial doubling time, days:", units = "user", cex = ct$caption$cex / cex, font = ct$caption$font, family = ct$caption$family)
      wu_lg_e <- strwidth("estimate", units = "user", cex = ct$estimate$cex / cex, font = ct$estimate$font, family = ct$caption$family)

      x_lg_cap <- usr[1L] + adj * (usr[2L] - usr[1L] - wu_lg_cap)
      x_lg_e <- x_lg_cap + adj * (wu_lg_cap - wu_lg_e) + 0.5 * wu_lg_e

      ## Estimates
      if (!is.null(ct$estimate)) {
        args <- list(
          text = c(sprintf("%.1f", elu[[1L]]), if (show_caption) "estimate"),
          side = 3,
          line = c(rep_len(line_e, n), if (show_caption) line_lg_e),
          at = c((t1 + t2) / 2, if (show_caption) x_lg_e),
          adj = 0.5,
          padj = 0
        )
        do.call(mtext, c(args, ct$estimate))
      }

      ## Confidence intervals
      if (show_tdoubling == 2L && !is.null(ct$ci)) {
        args <- list(
          text = c(sprintf("(%.1f, %.1f)", elu[[2L]], elu[[3L]]),
                   if (show_caption) sprintf("(%.3g%% CI)", 100 * level)),
          side = 3,
          line = c(rep_len(line_ci, n), if (show_caption) line_lg_ci),
          at = c((t1 + t2) / 2, if (show_caption) x_lg_e),
          adj = 0.5,
          padj = 0
        )
        do.call(mtext, c(args, ct$ci))
      }

      ## Caption
      if (show_caption) {
        args <- list(
          text = "initial doubling time, days:",
          side = 3,
          line = line_lg_cap,
          padj = 0
        )
        do.call(mtext, c(args, ct$caption))
      }
    }

    ## Plot (sub)title
    if (!is.null(ct1 <- control$title_main)) {
      line <- 0.5
      if (show_tdoubling > 0L) {
        line <- line_e + hin_to_lines(hin_e) + line
      }
      if (!is.null(ct2 <- control$title_sub)) {
        names(ct2) <- base::sub("\\.sub$", ".main", names(ct2))
        args <- list(main = sub[k], line = line)
        do.call(title, c(args, ct2))
        hin_sub <- strheight(sub[k], units = "inches", cex = ct2$cex.main, font = ct2$font.main, family = ct2$family)
        line <- line + hin_to_lines(hin_sub) + 0.25
      }
      args <- list(main = main[k], line = line)
      do.call(title, c(args, ct1))
    }

    ## Legend
    if (show_legend) {
      pty <- c("basic", "short", "long")
      ct <- control[sprintf("points_%s", pty)]

      f <- function(x) if (is.null(x)) NA else x
      ulf <- function(l) unlist(lapply(l, f), FALSE, FALSE)
      get_el <- function(el) ulf(lapply(ct, `[[`, el))

      args <- list(
        x = usr[2L] + 0.02 * (usr[2L] - usr[1L]),
        y = yu_to_y(usr[4L] - 0.02 * (usr[4L] - usr[3L])),
        xpd = NA,
        bty = "n",
        cex = 0.7,
        seg.len = 1
      )
      more_args <- list(
        legend = parse(
          text = sprintf("'%s,' ~ Delta * 't %s %g day%s'",
            c("obs", "obs", "obs", "pred"),
            c("=", "<", ">", "="),
            dt,
            if (dt != 1) "s" else ""
          )
        ),
        pch = c(get_el("pch"), NA),
        pt.bg = c(get_el("bg"), NA),
        lty = c(NA, NA, NA, f(control$lines$lty)),
        lwd = c(NA, NA, NA, f(control$lines$lwd)),
        col = c(get_el("col"), f(control$lines$col))
      )

      show_legend_item <- c(pty %in% levels(data$pty), TRUE)
      more_args <- lapply(more_args, `[`, show_legend_item)
      do.call(legend, c(args, more_args))
    }
  }

  invisible(NULL)
}

#' @import graphics
#' @importFrom grDevices colorRamp rgb
do_heat_plot <- function(cache, time_as, dt, origin, log,
                         panels_per_plot, inter_panel_space, control,
                         xlim, main, sub, xlab, ylab, ylab_outer, plab) {
  ### Set up ===================================================================

  range_time <- range(cache$time)
  shift <- range_time[1L]
  cache$time <- cache$time - shift

  range_log_r <- range(cache$estimate, na.rm = TRUE)
  diff_range_log_r <- diff(range_log_r)
  range_r <- exp(range_log_r)
  range_tdoubling <- log(2) / range_r[2:1]

  cache_split <- split(cache, cache$ts)
  N <- length(cache_split) # number of time series

  ## Axis limits (x)
  if (is.null(xlim)) {
    xlim <- range_time - shift
  } else if (inherits(xlim, "Date")) {
    xlim <- julian(xlim, origin = origin + shift)
  }

  ## Axis limits (y)
  ylim <- c(0, 1)

  ## Plot title
  if (is.null(main)) {
    main <- "Per capita growth rate, by time series"
  }

  ## Axis title (x)
  if (is.null(xlab)) {
    xlab <- switch(time_as,
      Date = "",
      sprintf("number of days since %s", origin + shift)
    )
  }

  ## Axis title (y)
  if (is.null(ylab)) {
    ylab <- "growth\nrate,\nper day"
  }
  if (is.null(ylab_outer)) {
    ylab_outer <- "doubling\ntime,\ndays"
  }

  ## Panel titles
  if (is.null(plab)) {
    plab <- names(cache_split)
  }

  ## Colour palette
  pal <- do.call(colorRamp, control$colorRamp)

  ## Utilities
  if (log) {
    to_unit <- function(x) (x - range_log_r[1L]) / diff_range_log_r
  } else {
    to_unit <- function(x) exp(x) / range_r[2L]
  }
  hin_to_lines <- function(hin) {
    diff(grconvertY(c(0, hin), "inches", "lines"))
  }


  ### Loop over plots ==========================================================

  ## Device layout
  L <- c(seq_len(panels_per_plot), rep_len(panels_per_plot + 1L, panels_per_plot))
  dim(L) <- c(panels_per_plot, 2L)
  layout(L, widths = c(par("din")[1L] - 2, 2))
  ## FIXME: `din` isn't updated until `plot.new` is called.
  ## Reliance on `din` here could cause unexpected behaviour in RStudio,
  ## e.g., if user resizes plot pane before `plot.new` call?

  ## Graphical parameters
  op <- par(oma = c(3.5, 0, 3.5, 0), xaxs = "i", yaxs = "i", cex = 1)
  on.exit(par(op))

  K <- 0L
  while (K < N) {
    ## Graphical parameters for heat map panels
    par(mar = c(0.5 * inter_panel_space, 2, 0.5 * inter_panel_space, 1))


    ### Loop over panels -------------------------------------------------------

    for (k in K + seq_len(min(panels_per_plot, N - K))) {
      d <- cache_split[[k]]

      plot.new()
      plot.window(xlim = xlim, ylim = ylim)
      usr <- par("usr")
      pin <- par("pin")

      ## Plot (sub)title
      if (k == K + 1L && !is.null(ct1 <- control$title_main)) {
        line <- 0.5
        if (!is.null(ct2 <- control$title_sub)) {
          names(ct2) <- base::sub("\\.sub$", ".main", names(ct2))
          args <- list(main = sub, line = line, xpd = NA)
          do.call(title, c(args, ct2))
          hin_sub <- strheight(sub, units = "inches", cex = ct2$cex.main, font = ct2$font.main, family = ct2$family)
          line <- line + hin_to_lines(hin_sub) + 0.25
        }
        args <- list(main = main, line = line, xpd = NA)
        do.call(title, c(args, ct1))
      }

      ## Background
      if (!is.null(ct <- control$rect_bg_panel)) {
        args <- list(
          xleft   = usr[1L],
          xright  = usr[2L],
          ybottom = usr[3L],
          ytop    = usr[4L]
        )
        do.call(rect, c(args, ct))
      }

      ## Pixels
      for (i in seq_len(nrow(d))) {
        rect(
          xleft   = d$time[i] - 0.5,
          xright  = d$time[i] + 0.5,
          ybottom = 0,
          ytop    = 1,
          border  = NA,
          col = rgb(pal(to_unit(d$estimate[i])), maxColorValue = 255)
        )
      }

      ## Panel title
      if (!is.null(ct1 <- control$title_plab)) {
        names(ct1) <- base::sub("\\.lab$", "", names(ct1))

        ## Pad from top left corner
        px <- 0.075 * (pin[2L] / pin[1L]) * (usr[2L] - usr[1L])
        py <- 0.075 * (usr[4L] - usr[3L])

        ## Underlay
        if (!is.null(ct2 <- control$rect_bg_plab)) {
          wu_plab <- strwidth(plab[k], units = "user", cex = ct1$cex, font = ct1$font, family = ct1$family)
          hu_plab <- strheight(plab[k], units = "user", cex = ct1$cex, font = ct1$font, family = ct1$family)
          args <- list(
            xleft   = usr[1L],
            xright  = usr[1L] + wu_plab + 2 * px,
            ybottom = usr[4L] - hu_plab - 2 * py,
            ytop    = usr[4L]
          )
          do.call(rect, c(args, ct2))
        }

        ## Text
        args <- list(
          x = usr[1L] + px,
          y = usr[4L] - py,
          labels = plab[k],
          adj = c(0, 1)
        )
        do.call(text, c(args, ct1))
      }
    }

    if (time_as == "Date") {
      ## Axis (x)
      Daxis(
        origin = origin + shift + 1,
        minor = control$axis_x_Date_minor,
        major = control$axis_x_Date_major,
        show_minor = !is.null(control$axis_x_Date_minor),
        show_major = !is.null(control$axis_x_Date_major)
      )
    } else {
      ## Axis (x)
      if (!is.null(ct <- control$axis_x_numeric)) {
        args <- list(side = 1)
        do.call(baxis, c(args, ct))
      }

      ## Axis title (x)
      if (!is.null(ct <- control$title_xlab)) {
        args <- list(xlab = xlab, line = 2.5)
        do.call(title, c(args, ct))
      }
    }

    ## Skip empty panels to get to the last
    while (k %% panels_per_plot > 0L) {
      plot.new()
      k <- k + 1L
    }

    ## Color scale
    par(mar = c(par("mar")[c(1L, 4L, 3L)], 8))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
    dy <- 0.001
    for (y in seq.int(0, 1, by = 2 * dy)) {
      rect(
        xleft   = 0,
        xright  = 1,
        ybottom = y - dy,
        ytop    = y + dy,
        border = NA,
        col = rgb(pal(y), maxColorValue = 255)
      )
    }

    if (!is.null(ct1 <- control$axis_y)) {
      par(new = TRUE, ylog = FALSE)
      if (log) {
        plot.window(xlim = c(0, 1), ylim = range_r, log = "y")
      } else {
        plot.window(xlim = c(0, 1), ylim = c(0, range_r[2L]), log = "")
      }

      ## Axis (y, inner)
      args <- list(side = 4, las = 1)
      do.call(baxis, c(args, ct1))

      par(new = TRUE, ylog = FALSE)
      plot.window(xlim = c(0, 1), ylim = range_tdoubling, log = "y")

      ## Axis (y, outer)
      args <- list(side = 4, las = 1)
      ct1$mgp <- ct1$mgp + 4
      do.call(baxis, c(args, ct1))

      if (!is.null(ct2 <- control$title_ylab)) {
        par(new = TRUE, ylog = FALSE)
        plot.window(xlim = c(0, 1), ylim = c(0, 1), log = "")
        names(ct2) <- sub("\\.lab$", "", names(ct2))
        ct2$adj <- NULL

        ## Axis title (y, inner)
        args <- list(
          x = 0,
          y = 1 + diff(grconvertY(c(0, 0.5), "lines", "user")),
          labels = ylab,
          adj = c(0, 0),
          xpd = NA
        )
        do.call(text, c(args, ct2))

        ## Axis title (y, outer)
        args$x <- diff(grconvertX(c(0, 5), "lines", "user"))
        args$labels <- ylab_outer
        do.call(text, c(args, ct2))
      }
    }

    K <- K + panels_per_plot
  }

  invisible(NULL)
}

#' Define plot options
#'
#' Sets parameters controlling the graphical output of \code{\link{plot.egf}}.
#' Here, \code{x}, \code{type}, \code{time_as}, and \code{dt} refer to the
#' so-named arguments of \code{\link{plot.egf}}.
#'
#' @param points
#'   A named \link{list} of arguments to \code{\link{points}},
#'   affecting the appearance of observed data.
#'   [\code{type != "rt2"} only.]
#' @param points_short,points_long
#'   Alternatives to \code{points} used for counts over intervals
#'   shorter or longer than \code{dt} days.
#'   [\code{type = "interval"} only.]
#' @param lines
#'   A named \link{list} of arguments to \code{\link{lines}},
#'   affecting the appearance of predicted curves.
#'   [\code{type != "rt2"} only.]
#' @param polygon
#'   A named \link{list} of arguments to \code{\link{polygon}},
#'   affecting the appearance of confidence bands on predicted curves.
#'   [\code{type != "rt2"} only.]
#' @param rect
#'   A named \link{list} of arguments to \code{\link{rect}},
#'   affecting the appearance of fitting windows.
#'   [\code{type != "rt2"} only.]
#' @param rect_bg_panel,rect_bg_plab
#'   A named \link{list} of arguments to \code{\link{rect}},
#'   affecting the appearance of panel backgrounds and panel title underlays.
#'   [\code{type = "rt2"} only.]
#' @param abline
#'   A named \link{list} of arguments to \code{\link{abline}},
#'   affecting the appearance of the line drawn at \code{y = 0}.
#'   [\code{type = "rt1"} only.]
#' @param segments
#'   A named \link{list} of arguments to \code{\link{segments}},
#'   affecting the appearance of line segments drawn at \code{y = r}
#'   when \code{x$model$curve} is \code{"exponential"}, \code{"logistic"},
#'   or \code{"richards"}.
#'   [\code{type = "rt1"} only.]
#' @param axis_x_Date_minor,axis_x_Date_major,axis_x_numeric,axis_y
#'   Named \link{list}s of arguments to \code{\link{axis}},
#'   affecting the appearance of plot axes. \code{axis_x_Date_*}
#'   are used for \code{time_as = "Date"}. \code{axis_x_numeric}
#'   is used for \code{time_as = "numeric"}.
#' @param box
#'   A named \link{list} of arguments to \code{\link{box}},
#'   affecting the appearance of the box drawn around the plot region.
#'   [\code{type != "rt2"} only.]
#' @param title_main,title_sub,title_xlab,title_ylab,title_plab
#'   Named \link{list}s of arguments to \code{\link{title}},
#'   affecting the appearance of plot, axis, and panel titles.
#' @param mtext_tdoubling_caption,mtext_tdoubling_estimate,mtext_tdoubling_ci
#'   Named \link{list}s of arguments to \code{\link{mtext}},
#'   affecting the appearance of initial doubling times printed
#'   in the top margin.
#'   [\code{type != "rt2"} only.]
#' @param colorRamp
#'   A named \link{list} of arguments to \code{\link{colorRamp}},
#'   defining heat map colour palette.
#'   [\code{type = "rt2"} only.]
#'
#' @details
#' Unsupported and unmodifiable options are silently discarded.
#' Modifiable options that are unspecified are assigned a default
#' value defined internally. \code{egf_plot_control()} returns
#' a \link{list} of default values for all modifiable options.
#'
#' Setting an argument to \code{\link{NULL}} has the effect of
#' suppressing the corresponding plot element. For example,
#' to suppress observed data, set \code{points} to \code{NULL},
#' and perhaps also \code{points_short} and \code{points_long}.
#'
#' @return
#' A named list.
#'
#' @export
egf_plot_control <- function(points, points_short, points_long,
                             lines, polygon, rect,
                             rect_bg_panel, rect_bg_plab,
                             abline, segments,
                             axis_x_Date_minor, axis_x_Date_major,
                             axis_x_numeric, axis_y, box,
                             title_main, title_sub,
                             title_xlab, title_ylab, title_plab,
                             mtext_tdoubling_caption,
                             mtext_tdoubling_estimate, mtext_tdoubling_ci,
                             colorRamp) {
  nf <- names(formals(egf_plot_control))
  values <- mget(nf, ifnotfound = NA, inherits = FALSE)
  defaults <- list(
    points = list(
      pch = 21,
      col = "#BBBBBB",
      bg = "#DDDDDD",
      cex = 1
    ),
    points_short = list(
      pch = 1,
      col = "#882255",
      bg = NA,
      cex = 1
    ),
    points_long = list(
      pch = 16,
      col = "#882255",
      bg = NA,
      cex = 1
    ),
    lines = list(
      lty = 1,
      lwd = 2.5,
      col = "#44AA99"
    ),
    polygon = list(
      col = "#44AA9960",
      border = NA,
      lty = 1,
      lwd = 1
    ),
    rect = list(
      col = "#DDCC7740",
      border = NA,
      lty = 1,
      lwd = 1
    ),
    rect_bg_panel = list(
      col = "black",
      border = NA,
      lty = 1,
      lwd = 1
    ),
    rect_bg_panel = list(
      col = "#00000080",
      border = NA,
      lty = 1,
      lwd = 1
    ),
    abline <- list(
      lty = 2,
      lwd = 1,
      col = "black"
    ),
    segments <- list(
      lty = 3,
      lwd = 2,
      col = "black"
    ),
    axis_x_Date_minor = list(
      mgp = c(3, 0.25, 0),
      lwd = 1,
      col = "black",
      lwd.ticks = 1,
      col.ticks = "black",
      tcl = -0.2,
      gap.axis = 0,
      col.axis = "black",
      cex.axis = 0.9,
      font.axis = 1,
      family = "sans",
      xpd = FALSE
    ),
    axis_x_Date_major = list(
      mgp = c(3, 1.25, 0),
      lwd = 1,
      col = "black",
      lwd.ticks = 1,
      col.ticks = "black",
      tcl = 0,
      gap.axis = 0,
      col.axis = "black",
      cex.axis = 1.2,
      font.axis = 1,
      family = "sans",
      xpd = FALSE
    ),
    axis_x_numeric = list(
      mgp = c(3, 0.7, 0),
      lwd = 1,
      col = "black",
      lwd.ticks = 1,
      col.ticks = "black",
      tcl = -0.5,
      gap.axis = 0,
      col.axis = "black",
      cex.axis = 0.9,
      font.axis = 1,
      family = "sans",
      xpd = FALSE
    ),
    axis_y = list(
      mgp = c(3, 0.7, 0),
      lwd = 1,
      col = "black",
      lwd.ticks = 1,
      col.ticks = "black",
      tcl = -0.5,
      gap.axis = NA,
      col.axis = "black",
      cex.axis = 0.9,
      font.axis = 1,
      family = "sans",
      xpd = FALSE
    ),
    box = list(
      bty = "l",
      lty = 1,
      lwd = 1,
      col = "black"
    ),
    title_main = list(
      adj = 0,
      col.main = "black",
      cex.main = 1,
      font.main = 2,
      family = "sans"
    ),
    title_sub = list(
      adj = 0,
      col.sub = "black",
      cex.sub = 0.75,
      font.sub = 2,
      family = "sans"
    ),
    title_xlab = list(
      adj = 0.5,
      col.lab = "black",
      cex.lab = 1,
      font.lab = 1,
      family = "sans"
    ),
    title_ylab = list(
      adj = 0.5,
      col.lab = "black",
      cex.lab = 1,
      font.lab = 1,
      family = "sans"
    ),
    title_plab <- list(
      col.lab = "white",
      cex.lab = 1,
      font.lab = 2,
      family = "sans"
    ),
    mtext_tdoubling_caption = list(
      col = "black",
      cex = 0.7,
      font = 1,
      family = "sans",
      adj = 1
    ),
    mtext_tdoubling_estimate = list(
      col = "black",
      cex = 0.7,
      font = 2,
      family = "sans"
    ),
    mtext_tdoubling_ci = list(
      col = "black",
      cex = 0.7,
      font = 1,
      family = "sans"
    ),
    colorRamp <- list(
      colors = c(
        "#364B9A", "#4A7BB7", "#6EA6CD", "#98CAE1",
        "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366",
        "#F67E4B", "#DD3D2D", "#A50026"
      ),
      bias = 1,
      space = "rgb",
      interpolate = "linear"
    )
  )

  lmerge <- function(value, default) {
    if (is.null(value)) {
      return(NULL)
    }
    if (!is.list(value) || length(value) == 0L) {
      return(default)
    }
    m <- match(names(value), names(default), 0L)
    default[m] <- value[m > 0L]
    default
  }
  control <- Map(lmerge, value = values, default = defaults)
  class(control) <- c("egf_plot_control", "list")
  control
}
