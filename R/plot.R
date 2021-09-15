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
#'   per capita growth rate (\code{"rt"}), and
#'   per capita growth rate \emph{as heat map} (\code{"rt_heat"}).
#'   \code{"rt_heat"} displays \code{per_plot} time series in each plot.
#'   The rest display one time series in each plot.
#' @param time_as
#'   A \link{character} string indicating how numeric times are displayed
#'   on the bottom axis. The options are:
#'   as is (\code{"numeric"})
#'   and with a calendar (\code{"Date"}).
#'   In the latter case, numeric times are interpreted as numbers of days
#'   since \code{1970-01-01 00:00:00}.
#' @param log
#'   A \link{logical} flag. If \code{TRUE}, then the dependent variable
#'   is plotted on a logarithmic scale.
#' @param dt
#'   A positive number indicating an observation interval.
#'   Predicted curves are evaluated on grids with this spacing.
#'   When \code{type = "interval"}, counts observed over shorter
#'   or longer intervals \code{dt0} are scaled by a factor of
#'   \code{dt / dt0} so that their scale matches that of the curves.
#'   Scaled counts can be highlighted via \code{control}.
#'   If \code{x} specifies a model with day of week effects
#'   (\code{x$model$day_of_week > 0}), then setting \code{dt}
#'   has no effect as it is set to 1 internally.
#' @param per_plot
#'   A positive integer giving the number of time series displayed in
#'   each plot.
#'   [\code{type = "rt_heat"} only.]
#' @param subset
#'   An expression to be evaluated in the combined model frame
#'   (see \code{\link{egf_combine_frames}}). It must evaluate
#'   to a valid index vector for the rows of the data frame
#'   (see \code{\link{[.data.frame}}), and thus fitting windows.
#'   Only indexed fitting windows are highlighted in plots.
#'   Only time series with indexed fitting windows are plotted.
#'   The default (\code{\link{NULL}}) is to plot all time series
#'   and display all fitting windows in each time series.
#' @param order
#'   An expression to be evaluated in the combined model frame
#'   (see \code{\link{egf_combine_frames}}), typically a
#'   call to \code{\link{order}}, determining the order in which
#'   time series are plotted. It must evaluate to a permutation of
#'   \code{\link{seq_len}(\link{nrow}(combined))}.
#'   The default (\code{\link{NULL}}) is equivalent to
#'   \code{\link{seq_len}(\link{nrow}(combined))}.
#' @param cache
#'   An \code{"egf_plot_cache"} object returned by a previous evaluation
#'   of \code{plot(x)}. Fitted and predicted values and standard errors
#'   stored in \code{cache} are reused if possible to avoid recomputation.
#' @param plot
#'   A \link{logical} flag. If \code{FALSE}, then nothing is plotted.
#'   Useful when only the returned \code{"egf_plot_cache"} object is desired.
#' @param show_predict
#'   An integer flag: 2 is to draw predicted curves with confidence bands,
#'   1 is draw predicted curves only, 0 is to draw neither.
#'   \link[=logical]{Logical} values are coerced to integer.
#'   [\code{type != "rt_heat"} only.]
#' @param show_tdoubling
#'   An integer flag: 2 is to print initial doubling time estimates in
#'   the top margin with confidence intervals, 1 is to print estimates
#'   only, 0 is to print neither. \link[=logical]{Logical} values are
#'   coerced to integer. Supported only if \code{x$model$curve} is
#'   \code{"exponential"}, \code{"logistic"}, or \code{"richards"}.
#'   [\code{type != "rt_heat"} only.]
#' @param level
#'   A number in the interval (0,1). This is the confidence level used
#'   when \code{show_predict = 2} or \code{show_tdoubling = 2}.
#'   [\code{type != "rt_heat"} only.]
#' @param control
#'   An \code{"\link{egf_plot_control}"} object controlling the appearance
#'   of most plot elements.
#' @param xlim,ylim
#'   \link[=numeric]{Numeric} vectors of length 2 specifying axis limits,
#'   which are recycled for all plots. If \code{time_as = "Date"}, then
#'   \code{xlim} can instead be a \link{Date} vector or a \link{character}
#'   vector coercible to Date via \code{\link{as.Date}(xlim)}. \code{ylim}
#'   is unused by \code{type = "rt_heat"}.
#' @param main,sub,xlab,ylab,ylab_outer,plab
#'   \link[=character]{Character} strings or expressions used
#'   to generate plot (\code{main}, \code{sub}), axis (\code{xlab},
#'   \code{ylab}, \code{ylab_outer}), and panel (\code{plab}) titles.
#'   \code{main}, \code{xlab}, and \code{ylab} are supported for all
#'   values of \code{type}.
#'   \code{sub} is unused by \code{type = "rt_heat"}.
#'   \code{plab} is used by \code{type = "rt_heat"} only.
#'   \code{ylab_outer} is used by \code{type = "rt(_heat)?"} only.
#'   When \code{type != "rt_heat"}, \code{main} and \code{sub} are evaluated
#'   in the combined model frame (see \code{\link{egf_combine_frames}})
#'   in order to generate unique (sub)titles for each plot.
#'   When \code{type = "rt_heat"}, \code{plab} is evaluated similarly
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
#' an error is thrown during plotting. Hence the cache will be
#' available temporarily via \code{\link{.Last.value}}. To ensure
#' that the cache is available permanently, assign the result of
#' the call to \code{\link{plot}} to a name:
#' \code{cache <- plot(x, \dots)}.
#'
#' Caching functionality must be used with care, as mismatch between
#' \code{x} and \code{cache} will not be detected. Constructions such
#' as \code{plot(y, cache = plot(x, \dots), \dots)}, where \code{x}
#' and \code{y} are different objects, should not be expected to produce
#' correct results.
#'
#' See topic \code{\link{egf_eval}} for details on nonstandard evaluation
#' of \code{subset}, \code{order}, \code{main}, \code{sub}, and \code{plab}.
#'
#' @export
#' @importFrom stats fitted predict complete.cases
plot.egf <- function(x,
                     type = c("interval", "cumulative", "rt", "rt_heat"),
                     time_as = c("Date", "numeric"),
                     log = TRUE,
                     dt = 1,
                     per_plot = min(6L, nlevels(x$frame$ts)),
                     subset = NULL,
                     order = NULL,
                     cache = NULL,
                     plot = TRUE,
                     show_predict = TRUE,
                     show_tdoubling = TRUE,
                     level = 0.95,
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
  stop_if_not_true_false(plot)
  type <- match.arg(type)
  if (x$model$day_of_week > 0L) {
    dt <- 1
  } else {
    stop_if_not_number(dt, "positive")
  }
  if (type == "rt_heat") {
    show_predict <- 1L
    show_tdoubling <- 0L
  } else {
    stop_if_not_true_false(show_predict, allow_numeric = TRUE)
    show_predict <- min(2L, max(0L, as.integer(show_predict))) # coercion to '0:2'
    if (x$model$curve %in% c("exponential", "logistic", "richards")) {
      stop_if_not_true_false(show_tdoubling, allow_numeric = TRUE)
      show_tdoubling <- min(2L, max(0L, as.integer(show_tdoubling))) # coercion to '0:2'
    } else {
      show_tdoubling <- 0L
    }
    if (any(c(show_predict, show_tdoubling) == 2L)) {
      stop_if_not_number_in_interval(level, 0, 1, "()")
    }
  }

  combined <- egf_combine_frames(x)
  subset <- egf_eval_subset(substitute(subset), combined, parent.frame())
  if (length(subset) == 0L) {
    stop("'subset' indexes zero fitting windows, so there is nothing to plot.")
  }
  order <- egf_eval_order(substitute(order), combined, parent.frame())
  subset <- order[order %in% subset]
  frame_windows <- x$frame_windows[subset, , drop = FALSE]

  lts <- as.character(unique(frame_windows$ts))
  frame_windows$ts <- factor(frame_windows$ts, levels = lts)

  reorder <- do.call(base::order, unname(frame_windows[c("ts", "start")]))
  subset <- subset[reorder]
  frame_windows <- frame_windows[reorder, , drop = FALSE]

  lw <- as.character(frame_windows$window)
  frame_windows$window <- factor(frame_windows$window, levels = lw)

  frame <- x$frame
  frame$ts <- factor(frame$ts, levels = lts)
  frame$window <- factor(frame$window, levels = lw)
  frame <- frame[!is.na(frame$ts), , drop = FALSE]

  if (plot) {
    time_as <- match.arg(time_as)
    stop_if_not_true_false(log)
    if (type == "rt_heat") {
      stop_if_not_integer(per_plot, "positive")
    }
    stopifnot(inherits(control, "egf_plot_control"))
    if (!is.null(xlim)) {
      if (is.character(xlim)) {
        xlim <- try(as.Date(xlim))
      }
      if (inherits(xlim, "Date")) {
        xlim <- julian(xlim)
      }
      stop_if_not(
        is.numeric(xlim),
        length(xlim) == 2L,
        is.finite(xlim),
        xlim[1L] < xlim[2L],
        m = "Invalid 'xlim'."
      )
    }
    if (!is.null(ylim)) {
      stop_if_not(
        is.numeric(ylim),
        length(ylim) == 2L,
        is.finite(ylim),
        ylim[1L] < ylim[2L],
        m = "Invalid 'ylim'."
      )
    }

    subset1 <- subset[match(lts, frame_windows$ts, 0L)]
    if (type == "rt_heat") {
      plab <- egf_eval_label(substitute(plab), combined, parent.frame())[subset1]
    } else {
      main <- egf_eval_label(substitute(main), combined, parent.frame())[subset1]
      sub <- egf_eval_label(substitute(sub), combined, parent.frame())[subset1]
    }
  }

  ## If necessary, initialize cache
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
    stopifnot(inherits(cache, "egf_plot_cache"))
  }
  nc <- names(cache)

  ## If necessary, augment 'cache' with predicted values
  ## of whatever is being plotted and standard errors
  if (show_predict > 0L) {
    what <- sub("_.*$", "", type)
    ok <- cache$var == sprintf("log(%s)", what) & !(show_predict == 2L & is.na(cache$se))
    required <- setdiff(lw, cache$window[ok])
    if (length(required) > 0L) {
      m <- match(required, x$frame_windows$window, 0L)
      time_split <- Map(seq.int,
        from = x$frame_windows$start[m],
        to = x$frame_windows$end[m],
        by = dt
      )
      px <- predict(x,
        what = what,
        time = unlist(time_split, FALSE, FALSE),
        window = rep.int(x$frame_windows$window[m], lengths(time_split)),
        log = TRUE,
        se = (show_predict == 2L)
      )
      if (show_predict == 1L) {
        px$se <- NA_real_
      }
      cache <- rbind(cache, px[nc])
    }
  }

  ## If necessary, augment 'cache' with fitted values of 'log(r)'
  ## and standard errors
  if (show_tdoubling > 0L) {
    ok <- cache$var == "log(r)" & !(show_tdoubling == 2L & is.na(cache$se))
    required <- setdiff(lw, cache$window[ok])
    if (length(required) > 0L) {
      fx <- fitted(x,
        par = "log(r)",
        link = TRUE,
        se = (show_tdoubling == 2L),
        .subset = (x$frame_windows$window %in% required)
      )
      if (show_tdoubling == 1L) {
        fx$se <- NA_real_
      }
      fx$time <- NA_real_
      names(fx)[match("top", names(fx), 0L)] <- "var"
      cache <- rbind(cache, fx[nc])
    }
  }

  ## Clean up
  o <- do.call(base::order, unname(cache[-5L]))
  cache <- cache[o, , drop = FALSE]
  i <- !duplicated(cache[1:4])
  cache <- cache[i, , drop = FALSE]
  row.names(cache) <- NULL
  class(cache) <- c("egf_plot_cache", "data.frame")

  ## If not plotting, then return
  if (!plot) {
    return(invisible(cache))
  }

  ## If plotting, then create an instruction to return
  ## 'cache' if the low level plot function throws an error
  cache_bak <- cache
  on.exit({
    message("Augmented 'cache' returned despite error ...")
    return(invisible(cache_bak))
  })

  ## Extract only those rows of 'cache' needed by the low level plot function
  cache[1:3] <- Map(factor, cache[1:3],
    levels = list(c(sprintf("log(%s)", what), if (show_tdoubling > 0L) "log(r)" else NULL), lts, lw)
  )
  i <- complete.cases(cache[1:3])
  cache <- cache[i, , drop = FALSE]

  ## Compute confidence intervals
  cache[c("lower", "upper")] <- list(NA_real_)
  i <- (show_predict == 2L & cache$var == sprintf("log(%s)", what)) | (show_tdoubling == 2L & cache$var == "log(r)")
  if (any(i)) {
    cache[i, c("lower", "upper")] <- do.call(do_wald, c(cache[i, c("estimate", "se"), drop = FALSE], list(level = level)))
  }

  ## Approximation of per capita growth rate from cumulative counts
  ## is much less noisy but sensitive to (unknown) initial value,
  ## so retrieve the predicted initial value and hope that it's reasonable
  if (type == "rt") {
    px <- predict(x,
      what = "cumulative",
      time = frame_windows$start,
      window = frame_windows$window,
      log = FALSE,
      se = FALSE
    )
    frame_windows$c0 <- px$estimate[match(frame_windows$window, px$window, 0L)]
  }

  if (type == "rt_heat") {
    # plot.egf.heat(
    #   cache = cache,
    #   time_as = time_as,
    #   log = log,
    #   dt = dt,
    #   per_plot = per_plot,
    #   control = control,
    #   xlim = xlim,
    #   main = main,
    #   sub = sub,
    #   xlab = xlab,
    #   ylab = ylab,
    #   ylab_outer = ylab_outer,
    #   plab = plab
    # )
  } else {
    plot.egf.curve(
      model = x$model,
      frame = frame,
      frame_windows = frame_windows,
      cache = cache,
      type = type,
      time_as = time_as,
      log = log,
      dt = dt,
      show_predict = show_predict,
      show_tdoubling = show_tdoubling,
      level = level,
      control = control,
      xlim = xlim,
      ylim = ylim,
      main = main,
      sub = sub,
      xlab = xlab,
      ylab = ylab,
      ylab_outer = ylab_outer
    )
  }

  ## Discard exit instructions if low level plot function runs without error
  on.exit()
  invisible(cache)
}

#' @import graphics
plot.egf.curve <- function(model, frame, frame_windows, cache,
                           type, time_as, log, dt,
                           show_predict, show_tdoubling,
                           level, control, xlim, ylim,
                           main, sub, xlab, ylab, ylab_outer) {
  ### Set up ===================================================================

  n <- nlevels(frame$ts)
  formula <- as.formula(call("~", as.name(type), quote(time)))
  xlim_bak <- xlim
  ylim_bak <- ylim
  xlab_bak <- xlab
  elu <- c("estimate", "lower", "upper")

  ## Plot titles
  if (is.null(main)) {
    if (show_predict == 0L) {
      main <- ""
    } else {
      if (model$curve %in% c("gompertz", "richards")) {
        substr(model$curve, 1L, 1L) <- toupper(substr(model$curve, 1L, 1L))
      }
      main <- paste("Fitted", model$curve, "model")
    }
    main <- rep_len(main, n)
  }

  ## Plot subtitles
  if (is.null(sub)) {
    sub <- levels(frame$ts)
  }

  ## Axis title (y)
  if (is.null(ylab)) {
    if (type == "rt") {
      ylab <- "per capita growth rate"
    } else {
      ylab <- paste(type, "incidence")
    }
  }
  if (is.null(ylab_outer)) {
    if (type == "rt") {
      ylab_outer <- "doubling time"
    } else {
      ylab_outer <- ""
    }
  }

  ## Utilities
  inches_to_lines <- function(inches, axis) {
    f <- switch(tolower(axis), x = grconvertX, y = grconvertY, stop("Oops"))
    diff(f(c(0, inches), "inches", "lines"))
  }
  unlog <- if (log) function(x) 10^x else identity

  ## Hack for 'points' loop
  control$points_basic <- control$points

  ## Graphical parameters
  op <- par(mar = c(3.5, 5 + 4 * (type == "rt" && log), 4, 1) + 0.1)
  on.exit(par(op))


  ### Loop over plots ==========================================================

  for (i in seq_len(n)) {
    ### Set up for plot --------------------------------------------------------

    data <- frame[unclass(frame$ts) == i, , drop = FALSE]
    data_windows <- frame_windows[unclass(frame_windows$ts) == i, , drop = FALSE]
    if (show_tdoubling > 0L) {
      cache_r <- cache[unclass(cache$ts) == i & cache$var == "log(r)", , drop = FALSE]
      cache_r[elu] <- exp(cache_r[elu])
    }
    if (show_predict > 0L) {
      cache_predict <- cache[unclass(cache$ts) == i & cache$var == sprintf("log(%s)", type), , drop = FALSE]
      cache_predict[elu] <- exp(cache_predict[elu])
    }

    data$dt <- c(NA, diff(data$time))
    data[[type]] <- switch(type,
      interval = c(NA, data$x[-1L]) * dt / data$dt,
      cumulative = c(0, cumsum(data$x[-1L])),
      rt = local({
        y <- rep_len(NA_real_, nrow(data))
        f <- function(x, x0) c(diff(log(cumsum(c(x0, x))), lag = 2) / 2, NA)
        split(y, data$window, drop = TRUE) <- Map(f,
          x = split(data$x, data$window, drop = TRUE),
          x0 = data_windows$c0
        )
        y[!is.finite(y)] <- NA
        if (log) {
          y[y == 0] <- NA
        }
        y
      })
    )
    if (type == "interval") {
      data$pty <- factor(sign(data$dt - dt), levels = c(0, -1, 1), labels = c("basic", "short", "long"))
    } else {
      data$pty <- factor("basic")
    }

    ## Axis limits (x)
    xlim <- xlim_bak
    if (is.null(xlim)) {
      xlim <- c(0, max(data$time) * 1.01)
    }

    ## Axis limits (y)
    ylim <- ylim_bak
    if (is.null(ylim)) {
      y <- data[[type]]
      if (type == "rt" && show_predict > 0L) {
        y <- c(y, unlist(cache_predict[elu], FALSE, FALSE))
      }
      y <- y[!is.na(y)]
      if (length(y) == 0L || all(y == 0)) {
        if (log) {
          ylim <- c(0.1, 1)
        } else {
          ylim <- c(0, 1)
        }
      } else {
        if (log) {
          if (type == "rt") {
            ry <- range(y)
            ylim <- exp(log(ry) + c(-1, 1) * 0.04 * diff(log(ry)))
            ylim[!is.finite(ylim)] <- ry[!is.finite(ylim)]
          } else {
            if (all(y == 0 | y >= 1) && any(y > 1)) {
              ylim <- exp(c(-0.04, 1.04) * log(max(y)))
            } else {
              ry <- range(y[y > 0])
              ylim <- exp(log(ry) + c(-1, 1) * 0.04 * diff(log(ry)))
              ylim[!is.finite(ylim)] <- ry[!is.finite(ylim)]
            }
          }
        } else {
          ylim <- c(0, 1.04 * max(y))
        }
      }
      if (log) {
        data[[type]][data[[type]] == 0] <- ylim[1L]
      }
    } else {
      if (log) {
        data[[type]][data[[type]] == 0] <- NA
      }
    }

    ## Axis title (x)
    xlab <- xlab_bak
    if (is.null(xlab)) {
      xlab <- switch(time_as, Date = "", numeric = "time")
    }


    ### Plot -------------------------------------------------------------------

    plot.new()
    plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i", log = if (log) "y" else "")

    ## Fitting windows
    if (!is.null(ct <- control$rect)) {
      args <- list(
        xleft   = data_windows$start,
        ybottom = ylim[1L],
        xright  = data_windows$end,
        ytop    = ylim[2L]
      )
      do.call(rect, c(args, ct))
    }

    ## Observed data
    for (s in levels(data$pty)) {
      if (!is.null(ct <- control[[paste0("points_", s)]])) {
        args <- list(formula = formula, data = data, subset = (data$pty == s))
        do.call(points, c(args, ct))
      }
    }

    ## Confidence bands on predicted curves
    if (show_predict == 2L && !is.null(ct <- control$polygon)) {
      for (px in split(cache_predict, cache_predict$window, drop = TRUE)) {
        if (type == "cumulative") {
          c0 <- data[[type]][match(px$time[1L], data$time, 0L)]
          px$lower <- c0 + px$lower - px$lower[1L]
          px$upper <- c0 + px$upper - px$upper[1L]
        }
        args <- list(
          x = c(px$time, rev(px$time)),
          y = c(px$lower, rev(px$upper))
        )
        do.call(polygon, c(args, ct))
      }
    }

    ## Predicted curves
    if (show_predict > 0L && !is.null(ct <- control$lines)) {
      for (px in split(cache_predict, cache_predict$window, drop = TRUE)) {
        if (type == "cumulative") {
          c0 <- data[[type]][match(px$time[1L], data$time, 0L)]
          px$estimate <- c0 + px$estimate - px$estimate[1L]
        }
        args <- list(formula = estimate ~ time, data = px)
        do.call(lines, c(args, ct))
      }
    }

    ## Asymptote
    if (type == "rt" && !is.null(ct <- control$segments)) {
      args <- list(
        x0 = data_windows$start,
        x1 = data_windows$end,
        y0 = cache_r$estimate,
        y1 = cache_r$estimate
      )
      do.call(segments, c(args, ct))
    }

    ## Box
    if (!is.null(ct <- control$box)) {
      do.call(box, ct)
    }

    if (time_as == "Date") {
      ## Axis (x)
      Daxis(
        side = 1,
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
        args <- list(xlab = xlab, line = 3)
        do.call(title, c(args, ct))
      }
    }

    if (!is.null(ct1 <- control$axis_y)) {
      line <- ct1$mgp[2L]

      ## Axis (y), inner
      args1 <- list(side = 2, las = 1, at = axTicks(side = 2))
      if (type != "rt" && max(args1$at) >= 1e+05) {
        args1$labels <- get_scientific_labels(args1$at)
      } else {
        args1$labels <- as.character(args1$at)
      }
      cex.axis <- get_fill_cex(args1$labels,
        target = max(0, 3.5 - ct1$mgp[2L]),
        units = "lines",
        font = ct1$font.axis,
        family = ct1$family
      )
      ct1$cex.axis <- min(cex.axis, ct1$cex.axis)
      do.call(baxis, c(args1, ct1))
      width <- max(strwidth(args1$labels, units = "inches", cex = ct1$cex.axis, font = ct1$font.axis, family = ct1$family))
      line <- line + inches_to_lines(width, "x") + 1

      ## Axis title (y), inner
      if (!is.null(ct2 <- control$title_ylab)) {
        args2 <- list(ylab = ylab, line = line)
        do.call(title, c(args2, ct2))
        width <- strheight(ylab, units = "inches", cex = ct2$cex.lab, font = ct2$font.lab, family = ct2$family)
        line <- line + inches_to_lines(width, "x") + 1
      }

      if (type == "rt" && log) {
        par(new = TRUE)
        plot.window(xlim = xlim, ylim = base::log(2) / ylim, xaxs = "i", yaxs = "i", log = "y")

        ## Axis (y), outer
        args1 <- list(side = 2, las = 1, at = axTicks(side = 2))
        args1$labels <- as.character(args1$at)
        ct1$mgp <- line + ct1$mgp
        do.call(baxis, c(args1, ct1))
        width <- max(strwidth(args1$labels, units = "inches", cex = ct1$cex.axis, font = ct1$font.axis, family = ct1$family))
        line <- ct1$mgp[2L] + inches_to_lines(width, "x") + 1

        ## Axis title (y), outer
        if (!is.null(ct2 <- control$title_ylab)) {
          args2 <- list(ylab = ylab_outer, line = line)
          do.call(title, c(args2, ct2))
        }

        par(new = TRUE)
        plot.window(xlim = xlim, ylim = ylim, log = "y")
      }
    }

    ## Initial doubling times
    if (show_tdoubling > 0L) {
      tdoubling <- base::log(2) / cache_r[elu]
      names(tdoubling) <- elu[c(1L, 3L, 2L)]

      s <- c("caption", "estimate", "ci")
      ct <- control[paste0("mtext_tdoubling_", s)]
      names(ct) <- s

      height_ci <- strheight("", units = "inches", cex = ct$ci$cex, font = ct$ci$font, family = ct$ci$family)
      height_e <- strheight("", units = "inches", cex = ct$estimate$cex, font = ct$estimate$font, family = ct$estimate$family)

      line_ci     <- 0.25
      line_e      <- line_ci    + inches_to_lines(height_ci, "y") + 0.15
      line_lg_ci  <- line_e     + inches_to_lines(height_e,  "y") + 0.5
      line_lg_e   <- line_lg_ci + inches_to_lines(height_ci, "y") + 0.15
      line_lg_cap <- line_lg_e  + inches_to_lines(height_e,  "y") + 0.25

      lg_cap <- "initial doubling time:"
      lg_e <- "estimate"
      lg_ci <- sprintf("(%.3g%% CI)", 100 * level)

      adj <- ct$caption$adj

      width_lg_cap <- strwidth(lg_cap, units = "user", cex = ct$caption$cex / par("cex"), font = ct$caption$font, family = ct$caption$family)
      width_lg_e <- strwidth(lg_e, units = "user", cex = ct$estimate$cex / par("cex"), font = ct$estimate$font, family = ct$estimate$family)
      width_lg_ci <- strwidth(lg_ci, units = "user", cex = ct$ci$cex / par("cex"), font = ct$ci$font, family = ct$ci$family)
      width_lg_ci <- max(width_lg_e, width_lg_ci)

      x_lg_cap <- xlim[1L] + adj * max(0, xlim[2L] - xlim[1L] - width_lg_cap)
      x_lg_e <- x_lg_ci <- x_lg_cap + adj * (width_lg_cap - width_lg_ci) + 0.5 * width_lg_ci

      show_caption <- !is.null(ct$caption)

      ## Estimates
      if (!is.null(ct$estimate)) {
        args <- list(
          text = c(sprintf("%.1f", tdoubling[["estimate"]]), if (show_caption) lg_e),
          side = 3,
          line = c(rep_len(line_e, nrow(data_windows)), if (show_caption) line_lg_e),
          at = c((data_windows$start + data_windows$end) / 2, if (show_caption) x_lg_e),
          adj = 0.5,
          padj = 0
        )
        do.call(mtext, c(args, ct$estimate))
      }

      ## Confidence intervals
      if (show_tdoubling == 2L && !is.null(ct$ci)) {
        args <- list(
          text = c(sprintf("(%.1f, %.1f)", tdoubling[["lower"]], tdoubling[["upper"]]), if (show_caption) lg_ci),
          side = 3,
          line = c(rep_len(line_ci, nrow(data_windows)), if (show_caption) line_lg_ci),
          at = c((data_windows$start + data_windows$end) / 2, if (show_caption) x_lg_ci),
          adj = 0.5,
          padj = 0
        )
        do.call(mtext, c(args, ct$ci))
      }

      ## Caption
      if (show_caption) {
        args <- list(
          text = lg_cap,
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
        line <- line_e + inches_to_lines(height_e, "y") + line
      }
      if (!is.null(ct2 <- control$title_sub)) {
        names(ct2) <- base::sub("\\.sub$", ".main", names(ct2))
        args <- list(main = sub[i], line = line)
        do.call(title, c(args, ct2))
        height <- strheight(sub[i], units = "inches", cex = ct2$cex.main, font = ct2$font.main, family = ct2$family)
        line <- line + inches_to_lines(height, "y") + 0.25
      }
      args <- list(main = main[i], line = line)
      do.call(title, c(args, ct1))
    }
  }

  invisible(NULL)
}

# #' @import graphics
# #' @importFrom grDevices colorRamp rgb
# plot.egf.heat <- function(cache, time_as, log, dt,
#                           per_plot, control,
#                           xlim, main, sub, xlab, ylab, ylab_outer, plab) {
#   ### Set up ===================================================================
#
#   range_log_r <- range(cache$estimate, na.rm = TRUE)
#   diff_range_log_r <- diff(range_log_r)
#   range_r <- exp(range_log_r)
#   range_tdoubling <- log(2) / range_r[2:1]
#
#   cache_split <- split(cache, cache$ts)
#   K <- length(cache_split) # number of time series
#
#   ## Axis limits (x)
#   if (is.null(xlim)) {
#     xlim <- range(cache$time)
#   } else if (inherits(xlim, "Date")) {
#     xlim <- julian(xlim)
#   }
#
#   ## Axis limits (y)
#   ylim <- c(0, 1)
#
#   ## Plot title
#   if (is.null(main)) {
#     main <- "Per capita growth rate, by time series"
#   }
#
#   ## Axis title (x)
#   if (is.null(xlab)) {
#     xlab <- switch(time_as, Date = "", numeric = "time")
#   }
#
#   ## Axis title (y)
#   if (is.null(ylab)) {
#     ylab <- "growth\nrate,\nper day"
#   }
#   if (is.null(ylab_outer)) {
#     ylab_outer <- "doubling\ntime,\ndays"
#   }
#
#   ## Panel titles
#   if (is.null(plab)) {
#     plab <- names(cache_split)
#   }
#
#   ## Colour palette
#   pal <- do.call(colorRamp, control$colorRamp)
#
#   ## Utilities
#   if (log) {
#     to_unit <- function(x) (x - range_log_r[1L]) / diff_range_log_r
#   } else {
#     to_unit <- function(x) exp(x) / range_r[2L]
#   }
#   hin_to_lines <- function(hin) {
#     diff(grconvertY(c(0, hin), "inches", "lines"))
#   }
#
#
#   ### Loop over plots ==========================================================
#
#   ## Device layout
#   L <- c(seq_len(per_plot), rep_len(per_plot + 1L, per_plot))
#   dim(L) <- c(per_plot, 2L)
#   layout(L, widths = c(par("din")[1L] - 2, 2))
#   ## FIXME: 'din' isn't updated until 'plot.new' is called.
#   ## Reliance on 'din' here could cause unexpected behaviour in RStudio,
#   ## e.g., if user resizes plot pane before 'plot.new' call?
#
#   ## Graphical parameters
#   op <- par(oma = c(3.5, 0, 3.5, 0), xaxs = "i", yaxs = "i", cex = 1)
#   on.exit(par(op))
#
#   k <- 0L
#   while (k < K) {
#     ## Graphical parameters for heat map panels
#     par(mar = c(0.5 * inter_panel_space, 2, 0.5 * inter_panel_space, 1))
#
#
#     ### Loop over panels -------------------------------------------------------
#
#     for (j in k + seq_len(min(per_plot, K - k))) {
#       d <- cache_split[[j]]
#
#       plot.new()
#       plot.window(xlim = xlim, ylim = ylim)
#       usr <- par("usr")
#       pin <- par("pin")
#
#       ## Plot (sub)title
#       if (j == k + 1L && !is.null(ct1 <- control$title_main)) {
#         line <- 0.5
#         if (!is.null(ct2 <- control$title_sub)) {
#           names(ct2) <- base::sub("\\.sub$", ".main", names(ct2))
#           args <- list(main = sub, line = line, xpd = NA)
#           do.call(title, c(args, ct2))
#           hin_sub <- strheight(sub, units = "inches", cex = ct2$cex.main, font = ct2$font.main, family = ct2$family)
#           line <- line + hin_to_lines(hin_sub) + 0.25
#         }
#         args <- list(main = main, line = line, xpd = NA)
#         do.call(title, c(args, ct1))
#       }
#
#       ## Background
#       if (!is.null(ct <- control$rect_bg_panel)) {
#         args <- list(
#           xleft   = usr[1L],
#           xright  = usr[2L],
#           ybottom = usr[3L],
#           ytop    = usr[4L]
#         )
#         do.call(rect, c(args, ct))
#       }
#
#       ## Pixels
#       for (i in seq_len(nrow(d))) {
#         rect(
#           xleft   = d$time[i] - 0.5,
#           xright  = d$time[i] + 0.5,
#           ybottom = 0,
#           ytop    = 1,
#           border  = NA,
#           col = rgb(pal(to_unit(d$estimate[i])), maxColorValue = 255)
#         )
#       }
#
#       ## Panel title
#       if (!is.null(ct1 <- control$title_plab)) {
#         names(ct1) <- base::sub("\\.lab$", "", names(ct1))
#
#         ## Pad from top left corner
#         px <- 0.075 * (pin[2L] / pin[1L]) * (usr[2L] - usr[1L])
#         py <- 0.075 * (usr[4L] - usr[3L])
#
#         ## Underlay
#         if (!is.null(ct2 <- control$rect_bg_plab)) {
#           wu_plab <- strwidth(plab[j], units = "user", cex = ct1$cex, font = ct1$font, family = ct1$family)
#           hu_plab <- strheight(plab[j], units = "user", cex = ct1$cex, font = ct1$font, family = ct1$family)
#           args <- list(
#             xleft   = usr[1L],
#             xright  = usr[1L] + wu_plab + 2 * px,
#             ybottom = usr[4L] - hu_plab - 2 * py,
#             ytop    = usr[4L]
#           )
#           do.call(rect, c(args, ct2))
#         }
#
#         ## Text
#         args <- list(
#           x = usr[1L] + px,
#           y = usr[4L] - py,
#           labels = plab[j],
#           adj = c(0, 1)
#         )
#         do.call(text, c(args, ct1))
#       }
#     }
#
#     if (time_as == "Date") {
#       ## Axis (x)
#       Daxis(
#         minor = control$axis_x_Date_minor,
#         major = control$axis_x_Date_major,
#         show_minor = !is.null(control$axis_x_Date_minor),
#         show_major = !is.null(control$axis_x_Date_major)
#       )
#     } else {
#       ## Axis (x)
#       if (!is.null(ct <- control$axis_x_numeric)) {
#         args <- list(side = 1)
#         do.call(baxis, c(args, ct))
#       }
#
#       ## Axis title (x)
#       if (!is.null(ct <- control$title_xlab)) {
#         args <- list(xlab = xlab, line = 2.5)
#         do.call(title, c(args, ct))
#       }
#     }
#
#     ## Skip empty panels to get to the last
#     while (j %% per_plot > 0L) {
#       plot.new()
#       j <- j + 1L
#     }
#
#     ## Color scale
#     par(mar = c(par("mar")[c(1L, 4L, 3L)], 8))
#     plot.new()
#     plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
#     dy <- 0.001
#     for (y in seq.int(0, 1, by = 2 * dy)) {
#       rect(
#         xleft   = 0,
#         xright  = 1,
#         ybottom = y - dy,
#         ytop    = y + dy,
#         border = NA,
#         col = rgb(pal(y), maxColorValue = 255)
#       )
#     }
#
#     if (!is.null(ct1 <- control$axis_y)) {
#       par(new = TRUE, ylog = FALSE)
#       if (log) {
#         plot.window(xlim = c(0, 1), ylim = range_r, log = "y")
#       } else {
#         plot.window(xlim = c(0, 1), ylim = c(0, range_r[2L]), log = "")
#       }
#
#       ## Axis (y, inner)
#       args <- list(side = 4, las = 1)
#       do.call(baxis, c(args, ct1))
#
#       par(new = TRUE, ylog = FALSE)
#       plot.window(xlim = c(0, 1), ylim = range_tdoubling, log = "y")
#
#       ## Axis (y, outer)
#       args <- list(side = 4, las = 1)
#       ct1$mgp <- ct1$mgp + 4
#       do.call(baxis, c(args, ct1))
#
#       if (!is.null(ct2 <- control$title_ylab)) {
#         par(new = TRUE, ylog = FALSE)
#         plot.window(xlim = c(0, 1), ylim = c(0, 1), log = "")
#         names(ct2) <- sub("\\.lab$", "", names(ct2))
#         ct2$adj <- NULL
#
#         ## Axis title (y, inner)
#         args <- list(
#           x = 0,
#           y = 1 + diff(grconvertY(c(0, 0.5), "lines", "user")),
#           labels = ylab,
#           adj = c(0, 0),
#           xpd = NA
#         )
#         do.call(text, c(args, ct2))
#
#         ## Axis title (y, outer)
#         args$x <- diff(grconvertX(c(0, 5), "lines", "user"))
#         args$labels <- ylab_outer
#         do.call(text, c(args, ct2))
#       }
#     }
#
#     k <- k + per_plot
#   }
#
#   invisible(NULL)
# }

#' Define plot options
#'
#' Sets parameters controlling the graphical output of \code{\link{plot.egf}}.
#' Here, \code{x}, \code{type}, \code{time_as}, and \code{dt} refer to the
#' so-named arguments of \code{\link{plot.egf}}.
#'
#' @param points
#'   A named \link{list} of arguments to \code{\link{points}},
#'   affecting the appearance of observed data.
#'   [\code{type != "rt_heat"} only.]
#' @param points_short,points_long
#'   Alternatives to \code{points} used for counts over intervals
#'   shorter or longer than \code{dt} days.
#'   [\code{type = "interval"} only.]
#' @param lines
#'   A named \link{list} of arguments to \code{\link{lines}},
#'   affecting the appearance of predicted curves.
#'   [\code{type != "rt_heat"} only.]
#' @param polygon
#'   A named \link{list} of arguments to \code{\link{polygon}},
#'   affecting the appearance of confidence bands on predicted curves.
#'   [\code{type != "rt_heat"} only.]
#' @param rect
#'   A named \link{list} of arguments to \code{\link{rect}},
#'   affecting the appearance of fitting windows.
#'   [\code{type != "rt_heat"} only.]
#' @param rect_bg_panel,rect_bg_plab
#'   A named \link{list} of arguments to \code{\link{rect}},
#'   affecting the appearance of panel backgrounds and panel title underlays.
#'   [\code{type = "rt_heat"} only.]
#' @param abline
#'   A named \link{list} of arguments to \code{\link{abline}},
#'   affecting the appearance of the line drawn at \code{y = 0}.
#'   [\code{type = "rt"} only.]
#' @param segments
#'   A named \link{list} of arguments to \code{\link{segments}},
#'   affecting the appearance of line segments drawn at \code{y = r}
#'   when \code{x$model$curve} is \code{"exponential"}, \code{"logistic"},
#'   or \code{"richards"}.
#'   [\code{type = "rt"} only.]
#' @param axis_x_Date_minor,axis_x_Date_major,axis_x_numeric,axis_y
#'   Named \link{list}s of arguments to \code{\link{axis}},
#'   affecting the appearance of plot axes. \code{axis_x_Date_*}
#'   are used for \code{time_as = "Date"}. \code{axis_x_numeric}
#'   is used for \code{time_as = "numeric"}.
#' @param box
#'   A named \link{list} of arguments to \code{\link{box}},
#'   affecting the appearance of the box drawn around the plot region.
#'   [\code{type != "rt_heat"} only.]
#' @param title_main,title_sub,title_xlab,title_ylab,title_plab
#'   Named \link{list}s of arguments to \code{\link{title}},
#'   affecting the appearance of plot, axis, and panel titles.
#' @param mtext_tdoubling_caption,mtext_tdoubling_estimate,mtext_tdoubling_ci
#'   Named \link{list}s of arguments to \code{\link{mtext}},
#'   affecting the appearance of initial doubling times printed
#'   in the top margin.
#'   [\code{type != "rt_heat"} only.]
#' @param colorRamp
#'   A named \link{list} of arguments to \code{\link{colorRamp}},
#'   defining heat map colour palette.
#'   [\code{type = "rt_heat"} only.]
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
