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
#'   \code{type != "rt2"} displays one time series per plot.
#'   \code{type = "rt2"} displays \code{per_plot} time series per plot.
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
#'   of \code{plot.egf(x)}. Predicted values and standard errors stored
#'   in \code{cache} will be reused (rather than recomputed) if possible.
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
#' @param per_plot
#'   A positive integer. This is the number of panels (time series)
#'   displayed in one plot.
#'   [\code{type = "rt2"} only.]
#' @param ips (\code{type = "rt2"} only.)
#'   A non-negative number specifying the space between panels
#'   ("interpanel space") in multipanel plots as a number of margin lines.
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
#'   \link[=character]{Character} strings or expressions used to generate
#'   plot (\code{main}, \code{sub}), axis (\code{xlab}, \code{ylab},
#'   \code{ylab_outer}), and panel (\code{plab}) titles. \code{main},
#'   \code{xlab}, and \code{ylab} are supported for all values of \code{type}.
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
#' @param align_right_tdoubling_caption
#'   A \link{logical} flag, used if \code{show_tdoubling > 0}
#'   to determine whether the caption explaining doubling time
#'   annotation should be aligned left or right in the top margin.
#'   [\code{type != "rt2"} only.]
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A \link[=data.frame]{data.frame} inheriting from \link{class}
#' \code{"egf_plot_cache"}. If argument \code{cache} was supplied
#' in the function call, then this data frame is the result of
#' augmenting \code{cache} with new computations.
#'
#' This object is returned \emph{in spite of} errors thrown
#' during plotting, to avoid waste of computation time.
#'
#' @details
#' Caching functionality must be used with care, as mismatch between \code{x}
#' and \code{cache} will \emph{not} be detected. Constructions such as
#' \code{plot(x2, cache = plot(x1))} should \emph{not} be expected to produce
#' correct results.
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
                     per_plot = 6L,
                     ips = 0.25,
                     control = egf_plot_control(),
                     xlim = NULL,
                     ylim = NULL,
                     main = NULL,
                     sub = NULL,
                     xlab = NULL,
                     ylab = NULL,
                     ylab_outer = NULL,
                     plab = NULL,
                     align_right_tdoubling_caption = TRUE,
                     ...) {
  ## FIXME: `order` _actually_ orders `levels(window)`,
  ## not `levels(ts)`, leading to some indexing gore ...

  type <- match.arg(type)
  stop_if_not_true_false(do_plot)
  stop_if_not_number_in_interval(dt, 0, Inf, "()")
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
    stop_if_not(
      inherits(control, "egf_plot_control"),
      m = "`control` must inherit from class \"egf_plot_control\". See `?egf_plot_control`."
    )

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

    order <- order_to_index(substitute(order), data = combined, enclos = parent.frame())
    if (type == "rt2") {
      plab <- label_to_character(substitute(plab), data = combined, enclos = parent.frame())
    } else {
      main <- label_to_character(substitute(main), data = combined, enclos = parent.frame())
      sub <- label_to_character(substitute(sub), data = combined, enclos = parent.frame())
    }

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
      stop_if_not_integer(per_plot, kind = "positive")
      stop_if_not_number_in_interval(ips, 0, Inf, "[)")
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
          ## FIXME: improve `wrap` to handle cases like this
          wrap(
            "Missing values preventing calculation of cumulative incidence. ",
            "These time series will not be displayed:"
          ),
          "\n\n", paste0("  ", tsl_subset[!keep], collapse = "\n"),
          .call = FALSE
        )
        tsl_subset <- tsl_subset[keep]
        frame$ts <- factor(frame$ts, levels = tsl_subset)
        subset <- subset & (tsl %in% tsl_subset)
        wl_subset <- wl[subset]
        frame$window <- factor(frame$window, levels = wl_subset)
      }
      if (all(!keep)) {
        stop("Nothing left to do ...", call. = FALSE)
      }
    }

    if (show_tdoubling > 0L) {
      stop_if_not_true_false(align_right_tdoubling_caption)
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
  on.exit({
    message("Augmented `cache` returned despite error ...")
    return(invisible(cache))
  })
  cache0 <- cache

  ## Extract only those rows of `cache` needed by the low level plot function
  cache0[nc[1:3]] <- Map(factor,
    x = cache0[nc[1:3]],
    levels = list(c(what, if (show_tdoubling > 0L) "log(r)"), tsl_subset, wl_subset)
  )
  i <- complete.cases(cache0[nc[1:3]])
  cache0 <- cache0[i, , drop = FALSE]

  ## Compute confidence intervals
  i <- (show_predict == 2L & cache0$var == what) | (show_tdoubling == 2L & cache0$var == "log(r)")
  cache0[c("lower", "upper")] <- NA_real_
  cache0[i, c("lower", "upper")] <- do_wald(
    estimate = cache0$estimate[i],
    se = cache0$se[i],
    level = level
  )

  ## Subset and order plot/axis/panel titles
  ## so that they correspond elementwise to `tsl_subset`
  subset <- order[order %in% subset]
  m <- match(tsl_subset, tsl[subset], 0L)

  if (type == "rt2") {
    stop("Oops")
    # do_heat_plot(
    #   cache = cache0,
    #   origin = origin,
    #   time_as = time_as,
    #   origin = attr(x$frame, "origin"),
    #   log = log,
    #   per_plot = per_plot,
    #   ips = ips,
    #   control = control,
    #   xlim = xlim,
    #   main = main,
    #   xlab = xlab,
    #   ylab = ylab,
    #   ylab_outer = ylab_outer,
    #   plab = plab[subset][m]
    # )
  } else {
    do_curve_plot(
      frame = frame,
      cache = cache0,
      type = type,
      time_as = time_as,
      dt = dt,
      origin = attr(x$frame, "origin"),
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
      ylab_outer = ylab_outer,
      align_right_tdoubling_caption
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
                          main, sub, xlab, ylab, ylab_outer,
                          align_right_tdoubling_caption) {
  ### Set up ===================================================================

  N <- nlevels(frame$ts)
  frame_split <- split(frame, frame$ts)
  cache_split <- split(cache, cache$ts)
  inv_log10 <- if (log) function(x) 10^x else identity
  formula <- as.formula(call("~", as.name(type), quote(time)))
  what <- switch(type, interval = "log_int_inc", cumulative = "log_cum_inc", "log_rt")

  ## Plot title
  if (is.null(main)) {
    if (curve %in% c("gompertz", "richards")) {
      substr(curve, 1L, 1L) <- toupper(substr(curve, 1L, 1L))
    }
    main <- rep_len(sprintf("Fitted %s model", curve), N)
  }

  ## Plot subtitle
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

  ## Plot margins
  mar <- switch(type,
    rt1 = c(3.5, 4 + (is.null(ylim) || ylim[2L] > 0) * 4, 4, 1) + 0.1,
    c(3.5, 5, 4, 1 + 5.5 * show_legend) + 0.1
  )
  op <- par(mar = mar)
  on.exit(par(op))


  ### Loop over plots ==========================================================

  for (k in seq_len(N)) {
    ### Set up for plot --------------------------------------------------------

    frame <- droplevels(frame_split[[k]])
    shift <- min(frame$time)
    wl <- levels(frame$window)
    n <- length(wl)

    cache <- droplevels(cache_split[[k]])
    cache_predict <- cache[cache$var == what, , drop = FALSE]
    cache_predict$time <- cache_predict$time - shift
    cache_predict_split <- split(cache_predict, cache_predict$window)
    cache_log_r <- cache[cache$var == "log(r)", , drop = FALSE]

    data <- data.frame(
      time = frame$time - shift,
      dt = c(NA, diff(frame$time))
    )
    data[[type]] <- switch(type,
      interval = c(NA, frame$x[-1L]) * dt / data$dt,
      cumulative = cumsum(c(0, frame$x[-1L])),
      rt1 = c(NA, diff(base::log(frame$x))) / data$dt
    )
    data[[type]][!is.finite(data[[type]])] <- NA

    t12 <- c(tapply(data$time, frame$window, range, simplify = FALSE))
    t1 <- vapply(t12, `[`, 0, 1L)
    t2 <- vapply(t12, `[`, 0, 2L)

    ## Axis limits (x)
    if (is.null(xlim)) {
      xlim0 <- c(0, max(data$time) * 1.04)
    } else if (inherits(xlim, "Date")) {
      xlim0 <- julian(xlim, origin = origin + shift)
    } else {
      xlim0 <- xlim
    }

    ## Axis limits (y)
    if (is.null(ylim)) {
      if (type == "rt1") {
        ylim0 <- range(data[[type]], na.rm = TRUE)
        ylim0[1L] <- max(-base::log(2), ylim0[1L])
        ylim0[2L] <- min( base::log(2), ylim0[2L])
        ylim0 <- ylim0 + c(-1, 1) * 0.04 * (ylim0[2L] - ylim0[1L])
      } else {
        ymax <- max(data[[type]], na.rm = TRUE)
        if (log) {
          zero <- ymax^-0.04
          data[[type]][data[[type]] == 0] <- zero
        } else {
          zero <- 0
        }
        ylim0 <- c(zero, if (log) ymax^1.04 else ymax * 1.04)
      }
    } else {
      ylim0 <- ylim
    }

    ## Axis title (x)
    if (is.null(xlab)) {
      xlab0 <- switch(time_as,
        Date = "",
        sprintf("number of days since %s", origin + shift)
      )
    } else {
      xlab0 <- xlab
    }

    ## Point highlighting according to observation interval
    data$pty <- "basic"
    if (type == "interval") {
      data$pty[data$dt < dt] <- "short"
      data$pty[data$dt > dt] <- "long"
    }
    data$pty <- factor(data$pty)
    control$points_basic <- control$points # hack


    ### Plot -------------------------------------------------------------------

    plot.new()
    plot.window(xlim = xlim0, ylim = ylim0, xaxs = "i", yaxs = "i", log = if (log) "y" else "")
    gp <- par(c("cex", "csi", "cxy", "pin", "usr"))

    ## Fitting windows
    if (!is.null(ct <- control$rect)) {
      args <- list(
        xleft   = t1,
        ybottom = inv_log10(gp$usr[3L]),
        xright  = t2,
        ytop    = inv_log10(gp$usr[4L])
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
        args <- list(
          x = c(pd$time, rev(pd$time)),
          y = c(exp(pd$lower), rev(exp(pd$upper)))
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

    ## Initial doubling times
    if (show_tdoubling > 0L) {
      elu <- base::log(2) / exp(cache_log_r[c("estimate", "upper", "lower")])
      ct <- control[sprintf("text_tdoubling_%s", c("estimate", "ci", "caption"))]
      names(ct) <- sub("^text_tdoubling_", "", names(ct))
      show_caption <- !is.null(ct$caption)

      ## Much ado about choosing user coordinates when log scale is in effect
      h_ci <- strheight("", units = "user", cex = ct$ci$cex,       font = ct$ci$font)
      h_e  <- strheight("", units = "user", cex = ct$estimate$cex, font = ct$estimate$font)
      y_ci <- add_lines_to_user(0.25, inv_log10(gp$usr[4L]), log)
      y_e  <- add_lines_to_user(0.15, add_height_to_user(h_ci, y_ci, log), log)
      y_cap_ci <- add_lines_to_user(0.5,  add_height_to_user(h_e, y_e, log), log)
      y_cap_e  <- add_lines_to_user(0.15, add_height_to_user(h_ci, y_cap_ci, log), log)
      w_e <- strwidth("estimate", units = "user", cex = ct$estimate$cex, font = ct$estimate$font)
      if (align_right_tdoubling_caption) {
        x_cap <- gp$usr[2L] - 0.5 * w_e
      } else {
        x_cap <- gp$usr[1L] + 0.5 * w_e
      }

      if (!is.null(ct$estimate)) {
        args <- list(
          x = c((t1 + t2) / 2,
                if (show_caption) x_cap),
          y = c(rep_len(y_e, n),
                if (show_caption) y_cap_e),
          labels = c(sprintf("%.1f", elu[[1L]]),
                     if (show_caption) "estimate"),
          adj = c(0.5, 0),
          xpd = TRUE
        )
        do.call(text, c(args, ct$estimate))
      }

      if (show_tdoubling == 2L && !is.null(ct$ci)) {
        args <- list(
          x = c((t1 + t2) / 2,
                if (show_caption) x_cap),
          y = c(rep_len(y_ci, n),
                if (show_caption) y_cap_ci),
          labels = c(sprintf("(%.1f, %.1f)", elu[[2L]], elu[[3L]]),
                     if (show_caption) sprintf("(%.3g%% CI)", 100 * level)),
          adj = c(0.5, 0),
          xpd = TRUE
        )
        do.call(text, c(args, ct$ci))
      }

      if (show_caption) {
        args <- list(
          x = gp$usr[as.integer(align_right_tdoubling_caption) + 1L],
          y = add_lines_to_user(0.25, add_height_to_user(h_e, y_cap_e, log), log),
          labels = "initial doubling time, days:",
          adj = c(as.numeric(align_right_tdoubling_caption), 0),
          xpd = TRUE
        )
        do.call(text, c(args, ct$caption))
      }
    }

    ## Box
    if (!is.null(ct <- control$box)) {
      do.call(box, ct)
    }

    if (time_as == "Date") {
      ## Axis (x)
      daxis(
        origin = origin + shift + 1,
        minor = control$axis_x_Date_minor,
        major = control$axis_x_Date_major
      )
    } else {
      ## Axis (x)
      if (!is.null(ct1 <- control$axis_x_numeric)) {
        args <- list(side = 1)
        do.call(baxis, c(args, ct1))

        ## Axis title (x)
        if (!is.null(ct2 <- control$title_xlab)) {
          args <- list(xlab = xlab0, line = 2.5)
          do.call(title, c(args, ct2))
        }
      }
    }

    ## Axis (y)
    if (!is.null(ct1 <- control$axis_y)) {
      args <- list(
        side = 2,
        at = axTicks(side = 2),
        las = 1
      )
      if (type != "rt1" && max(args$at) >= 1e05) {
        args$labels <- get_yax_labels(args$at)
        ct1$cex.axis <- min(ct1$cex.axis, get_yax_cex(args$labels, mex = 3.5 - ct1$mgp[2L]))
      }
      do.call(baxis, c(args, ct1))
      if (type == "rt1" && gp$usr[4L] > 0) {
        tdoubling <- c(1:5, 10, 20, 50, 100)
        args <- list(
          side = 2,
          a = 0,
          b = gp$usr[4L],
          at = log(2) / tdoubling,
          labels = tdoubling,
          las = 1
        )
        ct1_outer <- ct1
        ct1_outer$mgp <- ct1$mgp + 4
        do.call(baxis, c(args, ct1_outer))
      }

      ## Axis title (y)
      if (!is.null(ct2 <- control$title_ylab)) {
        if (type == "rt1") {
          mlw <- max(strwidth(axTicks(side = 2),
            units = "inches",
            cex = ct1$cex.axis,
            font = ct1$font.axis
          ))
          line0 <- 0.5 + mlw / gp$csi + ct1$mgp[2L]
        } else {
          line0 <- 4
        }
        args <- list(ylab = ylab, line = line0)
        do.call(title, c(args, ct2))
        if (type == "rt1" && gp$usr[4L] > 0) {
          mlw <- max(strwidth(tdoubling,
            units = "inches",
            cex = ct1$cex.axis,
            font = ct1$font.axis
          ))
          tw <- strwidth(ylab_outer,
            units = "inches",
            cex = ct2$cex.lab,
            font = ct2$font.lab
          )
          ## `tw` relative to 0.8 times the distance from 0 to `usr[4]`
          rtw <- (tw * (gp$usr[4L] - gp$usr[3L]) / gp$pin[2L]) /
            (0.8 * (gp$usr[4L] - max(0, gp$usr[3L])))
          text(
            ## FIXME: Unexpected behavior
            x = gp$usr[1L] - gp$cxy[1L] * (2.5 + mlw / gp$csi + ct1_outer$mgp[2L]),
            y = mean(c(max(0, gp$usr[3L]), gp$usr[4L])),
            labels = ylab_outer,
            adj = c(0.5, 0),
            srt = 90,
            xpd = NA,
            col = ct2$col.lab,
            cex = ct2$cex.lab / max(1, rtw),
            font = ct2$font.lab
          )
        }
      }
    }

    ## Plot (sub)title
    if (!is.null(ct1 <- control$title_main)) {
      if (show_tdoubling > 0L) {
        line0 <- user_to_lines(add_lines_to_user(0.5, add_height_to_user(h_e, y_e, log), log), log)
      } else {
        line0 <- 0.5
      }
      if (!is.null(sub) && !is.null(ct2 <- control$title_sub)) {
        names(ct2) <- base::sub("\\.sub$", ".main", names(ct2))
        th <- strheight(sub[k], units = "inches", cex = ct2$cex.main, font = ct2$font.main)
        args <- list(main = sub[k], line = line0)
        do.call(title, c(args, ct2))
        line0 <- line0 + th / gp$csi + 0.25
      }
      args <- list(main = main[k], line = line0)
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
        x = gp$usr[2L] + 0.02 * (gp$usr[2L] - gp$usr[1L]),
        y = inv_log10(gp$usr[4L] - 0.02 * (gp$usr[4L] - gp$usr[3L]), log),
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
      do.call(legend, c(args, lapply(more_args, `[`, show_legend_item)))
    }
  }

  invisible(NULL)
}

# #' @import graphics
# #' @importFrom grDevices colorRamp rgb
# do_heat_plot <- function(cache, origin, time_as, log, per_plot, control,
#                          xlim, main, xlab, ylab, ylab_outer, plab) {
#   tr <- range(cache$time) # t range
#   shift <- tr[1L]
#   cache$time <- cache$time - shift
#
#   lrr <- range(cache$estimate, na.rm = TRUE) # log(r(t)) range
#   dlrr <- lrr[2L] - lrr[1L]
#   rr <- exp(lrr) # r(t) range
#
#   to_unit <- function(x, log) {
#     if (log) (x - lrr[1L]) / dlrr else exp(x) / rr[2L]
#   }
#   inv_log10 <- if (log) function(x) 10^x else identity
#
#   cache_split <- split(cache, cache$ts)
#   tsl <- levels(cache$ts)
#   N <- length(tsl)
#
#   ## Heat map color palette
#   pal <- do.call(colorRamp, control$colorRamp)
#   ips <- control$ips
#
#   ## Axis limits
#   if (is.null(xlim)) {
#     xlim <- tr - shift
#   } else {
#     if (inherits(xlim, "Date")) {
#       xlim <- julian(xlim, origin = origin + shift)
#     }
#   }
#   ylim <- c(0, 1)
#
#   ## Titles
#   if (is.null(main)) {
#     main <- "Per capita growth rate, by time series"
#   }
#   if (is.null(xlab)) {
#     xlab <- switch(time_as,
#       Date = "",
#       sprintf("number of days since %s", origin + shift)
#     )
#   }
#   if (is.null(ylab)) {
#     ylab <- "growth\nrate,\nper day"
#   }
#   if (is.null(ylab_outer)) {
#     ylab_outer <- "doubling\ntime,\ndays"
#   }
#   if (is.null(plab)) {
#     plab <- tsl
#   }
#
#   op <- par(oma = c(3.5, 0, 3.5, 0))
#   on.exit(par(op))
#
#   ### Loop over plots =====================================
#
#   K <- 0L
#   while (K < N) {
#     L <- c(seq_len(per_plot), rep_len(per_plot + 1L, per_plot))
#     dim(L) <- c(per_plot, 2L)
#     ## FIXME: Graphical parameter `din` isn't updated until
#     ## `plot.new()` is called or R session is restarted...
#     ## Use here could cause unexpected behavior in RStudio?
#     layout(L, widths = c(par("din")[1L] - 1.3, 1.3))
#     par(mar = c(0.5 * ips, 2, 0.5 * ips, 1))
#
#     ### Loop over panels ==================================
#
#     for (k in K + seq_len(min(per_plot, N - K))) {
#       d <- cache_split[[k]]
#       plot.new()
#       plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i")
#       gp <- par(c("cex", "csi", "cxy", "pin", "usr"))
#
#       ## Panel background
#       if (!is.null(control$rect$bg)) {
#         l <- list(
#           xleft   = gp$usr[1L],
#           xright  = gp$usr[2L],
#           ybottom = gp$usr[3L],
#           ytop    = gp$usr[4L]
#         )
#         do.call(rect, c(l, control$rect$bg))
#       }
#
#       ## Panel pixels
#       for (i in seq_len(nrow(d))) {
#         rect(
#           xleft   = d$time[i] - 0.5,
#           xright  = d$time[i] + 0.5,
#           ybottom = 0,
#           ytop    = 1,
#           border = NA,
#           col = rgb(pal(to_unit(d$estimate[i], log)), maxColorValue = 255)
#         )
#       }
#
#       ## Panel title
#       if (!is.null(ct <- control$title$plab)) {
#         names(ct) <- sub("\\.lab$", "", names(ct))
#         ct$cex <- gp$cex * ct$cex
#
#         ## Pad from top left corner
#         px <- (gp$usr[2L] - gp$usr[1L]) * (0.075 * gp$pin[2L] / gp$pin[1L])
#         py <- 0.075 * (gp$usr[4L] - gp$usr[3L])
#
#         ## Background
#         if (!is.null(control$rect$plab)) {
#           tw <- strwidth( plab[k], units = "user", cex = ct$cex, font = ct$font)
#           th <- strheight(plab[k], units = "user", cex = ct$cex, font = ct$font)
#           l <- list(
#             xleft   = gp$usr[1L],
#             xright  = gp$usr[1L] + tw + 2 * px,
#             ybottom = gp$usr[4L] - th - 2 * py,
#             ytop    = gp$usr[4L]
#           )
#           do.call(rect, c(l, control$rect$plab))
#         }
#
#         ## Text
#         l <- list(
#           x = gp$usr[1L] + px,
#           y = gp$usr[4L] - py,
#           labels = plab[k],
#           adj = c(0, 1)
#         )
#         do.call(text, c(l, ct))
#       }
#
#       ## Plot title
#       if (k == K + 1L && !is.null(ct <- control$title$main)) {
#         ct$cex.main <- ct$cex.main / gp$cex
#         l <- list(main = main, line = 0.5, xpd = NA)
#         do.call(title, c(l, ct))
#       }
#     }
#
#     ## Axis (x)
#     if (!is.null(ct <- control$axis$x)) {
#       if (time_as == "Date") {
#         ct$minor$cex.axis <- ct$minor$cex.axis / gp$cex
#         ct$major$cex.axis <- ct$major$cex.axis / gp$cex
#         daxis(
#           origin = origin + shift + 1,
#           minor = ct$minor,
#           major = ct$major
#         )
#       } else {
#         ct$cex.axis <- ct$cex.axis / gp$cex
#         l <- list(side = 1)
#         do.call(baxis, c(l, ct))
#       }
#     }
#
#     ## Axis title (x)
#     if (!is.null(ct <- control$title$xlab)) {
#       ct$cex.lab <- ct$cex.lab / gp$cex
#       l <- list(xlab = xlab, line = 2.5, xpd = NA)
#       do.call(title, c(l, ct))
#     }
#
#     ## Skip empty panels to get to the last
#     for (i in seq_len((per_plot - (k %% per_plot)) %% per_plot)) {
#       plot.new()
#     }
#
#     ## Color scale
#     par(mar = c(0.5 * ips, 1, 0.5 * ips, 8))
#     plot.new()
#     plot.window(
#       xlim = c(0, 1),
#       ylim = c(0, 1),
#       xaxs = "i",
#       yaxs = "i"
#     )
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
#     if (!is.null(ct <- control$axis$y)) {
#       par(new = TRUE)
#       plot.window(
#         xlim = c(0, 1),
#         ylim = if (log) rr else c(0, rr[2L]),
#         xaxs = "i",
#         yaxs = "i",
#         log = if (log) "y" else ""
#       )
#       gp <- par(c("cex", "csi", "cxy", "usr"))
#       ct$cex.axis <- ct$cex.axis / gp$cex
#
#       ## Axis (y, inner)
#       l <- list(side = 4, las = 1)
#       do.call(baxis, c(l, ct))
#
#       ## Axis (y, outer)
#       ct$mgp <- ct$mgp + c(4.5, 4, 4.5) # FIXME: Behaving unexpectedly...
#       tdoubling <- c(1:5, 10, 20, 50, 100)
#       l <- list(
#         side = 4,
#         at = base::log(2) / tdoubling,
#         labels = tdoubling,
#         las = 1
#       )
#       do.call(baxis, c(l, ct))
#
#       if (!is.null(ct <- control$title$y)) {
#         names(ct) <- sub("\\.lab$", "", names(ct))
#         y0 <- add_lines_to_user(0.5, inv_log10(gp$usr[4L]), log)
#
#         ## Axis title (y, inner)
#         l <- list(
#           x = gp$usr[2L],
#           y = y0,
#           labels = ylab,
#           adj = c(0, 0),
#           xpd = NA
#         )
#         do.call(text, c(l, ct))
#
#         ## Axis title (y, outer)
#         l <- list(
#           x = gp$usr[2L] + 4 * gp$cxy[1L],
#           y = y0,
#           labels = ylab_outer,
#           adj = c(0, 0),
#           xpd = NA
#         )
#         do.call(text, c(l, ct))
#       }
#     }
#     K <- K + per_plot
#   }
#
#   invisible(NULL)
# }
