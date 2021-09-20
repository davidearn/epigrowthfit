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
#'   In the latter case, horizontal user coordinates on measure time in days
#'   since \code{1970-01-01 00:00:00}.
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
#' @param log
#'   A \link{logical} flag. If \code{TRUE}, then the dependent variable
#'   is plotted on a logarithmic scale.
#' @param zero
#'   A positive number indicating a line on which to plot zeros
#'   when \code{log = TRUE}. \code{\link{NA}} is to place zeros
#'   on the bottom axis. \code{\link{NULL}} is to suppress zeros.
#'   [\code{type = "interval"} and \code{"cumulative"} only].
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
#' @param per_plot
#'   A positive integer giving the number of time series displayed
#'   in each plot.
#'   [\code{type = "rt_heat"} only.]
#' @param control
#'   An \code{"\link{egf_plot_control}"} object controlling the appearance
#'   of most plot elements.
#' @param cache
#'   An \code{"egf_plot_cache"} object returned by a previous evaluation
#'   of \code{plot(x)}. Fitted and predicted values and standard errors
#'   stored in \code{cache} are reused if possible to avoid recomputation.
#' @param plot
#'   A \link{logical} flag. If \code{FALSE}, then nothing is plotted.
#'   Useful when only the returned \code{"egf_plot_cache"} object is desired.
#' @param subset
#'   An expression to be evaluated in the combined model frame
#'   (see \code{\link{egf_combine_frames}}). It must evaluate
#'   to a valid index vector for the rows of the data frame
#'   (see \code{\link{[.data.frame}}), and thus fitting windows.
#'   Only time series with indexed fitting windows are plotted.
#'   Only indexed fitting windows are highlighted in plots.
#'   The default (\code{\link{NULL}}) is to plot all time series
#'   and display all fitting windows in each time series.
#' @param order
#'   An expression to be evaluated in the combined model frame
#'   (see \code{\link{egf_combine_frames}}), typically a
#'   call to \code{\link{order}}, determining the order in which
#'   time series are plotted. It must evaluate to a permutation
#'   of \code{\link{seq_len}(\link{nrow}(combined))}.
#'   The default (\code{\link{NULL}}) is the original ordering.
#' @param xlim,ylim
#'   \link[=numeric]{Numeric} vectors of length 2 specifying axis limits,
#'   which are recycled for all plots. If \code{time_as = "Date"},
#'   then \code{xlim} can instead be a \link{Date} vector or any vector
#'   coercible to Date via \code{\link{as.Date}(xlim)}.
#'   \code{ylim} is unused by \code{type = "rt_heat"}.
#' @param main,sub,xlab,ylab,ylab_outer,plab
#'   \link[=character]{Character} strings or expressions used
#'   to generate plot (\code{main}, \code{sub}), axis (\code{xlab},
#'   \code{ylab}, \code{ylab_outer}), and panel (\code{plab}) titles.
#'   \code{main}, \code{sub}, \code{xlab}, and \code{ylab} are supported
#'   for all values of \code{type}.
#'   \code{plab} is used by \code{type = "rt_heat"} only.
#'   \code{ylab_outer} is used by \code{type = "rt"} and \code{"rt_heat"} only.
#'   When \code{type != "rt_heat"}, \code{main} and \code{sub} are
#'   evaluated in the combined model frame
#'   (see \code{\link{egf_combine_frames}})
#'   in order to generate unique (sub)titles for each plot.
#'   When \code{type = "rt_heat"}, \code{plab} is evaluated
#'   similarly in order to generate unique titles for each panel.
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
#' an error is thrown during plotting. To ensure that the cache is
#' preserved, assign the result of the call to \code{\link{plot}}
#' to a name: \code{cache <- plot(x, \dots)}.
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
#' @examples
#' example("egf", "epigrowthfit")
#'
#' par(mar = c(3.5, 5, 5, 1))
#' control <- egf_plot_control(tdoubling = list(legend = list(cex = 0.8), estimate = list(cex = 0.8, font = 2), ci = list(cex = 0.8)))
#' plot(object, type = "interval", show_predict = 2L, show_tdoubling = 2L, control = control)
#' plot(object, type = "cumulative", main = "Fitted exponential model", sub = paste("Country", country))
#'
#' par(mar = c(3.5, 9.5, 5, 1))
#' plot(object, type = "rt", subset = country %in% LETTERS[4:6])
#'
#' par(mar = c(0.1, 1, 0.1, 1), oma = c(3.5, 0, 0.5, 9.5))
#' plot(object, type = "rt_heat", per_plot = 10L)
#'
#' @export
#' @importFrom stats fitted predict complete.cases
plot.egf <- function(x,
                     type = c("interval", "cumulative", "rt", "rt_heat"),
                     time_as = c("Date", "numeric"),
                     dt = 1,
                     log = TRUE,
                     zero = NA,
                     show_predict = TRUE,
                     show_tdoubling = FALSE,
                     level = 0.95,
                     per_plot = min(6L, nlevels(x$frame$ts)),
                     control = egf_plot_control(),
                     cache = NULL,
                     plot = TRUE,
                     subset = NULL,
                     order = NULL,
                     xlim = NULL,
                     ylim = NULL,
                     main = NULL,
                     sub = NULL,
                     xlab = NULL,
                     ylab = NULL,
                     ylab_outer = NULL,
                     plab = NULL,
                     ...) {
  ## Validate arguments --------------------------------------------------------

  type <- match.arg(type)
  if (x$model$day_of_week > 0L) {
    dt <- 1
  } else {
    stop_if_not_number(dt, "positive")
  }
  if (type == "rt_heat") {
    show_predict <- 1L
    show_tdoubling <- 0L
    show_asymptote <- 0L
  } else {
    stop_if_not_true_false(show_predict, allow_numeric = TRUE)
    show_predict <- min(2L, max(0L, as.integer(show_predict)))
    if (x$model$curve %in% c("exponential", "logistic", "richards")) {
      stop_if_not_true_false(show_tdoubling, allow_numeric = TRUE)
      show_tdoubling <- min(2L, max(0L, as.integer(show_tdoubling)))
      show_asymptote <- as.integer(x$model$curve != "exponential")
    } else {
      show_tdoubling <- 0L
      show_asymptote <- 0L
    }
  }
  stop_if_not_true_false(plot)
  if (plot) {
    time_as <- match.arg(time_as)
    stop_if_not_true_false(log)
    if (!is.null(zero)) {
      if (grepl("^rt(_heat)?$", type)) {
        zero <- NULL
      } else if (is.na(zero)) {
        zero <- as.numeric(zero)
      } else {
        stop_if_not_number(zero, "positive")
      }
    }
    if (any(c(show_predict, show_tdoubling) == 2L)) {
      stop_if_not_number_in_interval(level, 0, 1, "()")
    }
    if (type == "rt_heat") {
      stop_if_not_integer(per_plot, "positive")
    }
    stopifnot(inherits(control, "egf_plot_control"))
    if (!is.null(cache)) {
      stopifnot(inherits(cache, "egf_plot_cache"))
    }
    if (!is.null(xlim)) {
      if (!is.numeric(xlim)) {
        xlim <- try(julian(as.Date(xlim)))
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
  }


  ## Subset and order time series, fitting windows -----------------------------

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


  ## Finalize annotation -------------------------------------------------------

  if (plot) {
    n <- nlevels(frame$ts)
    subset1 <- subset[match(lts, frame_windows$ts, 0L)]
    if (type == "rt_heat") {
      nplot <- as.integer(ceiling(n / per_plot))
      if (is.null(main)) {
        main <- ""
      }
      if (is.null(sub)) {
        sub <- ""
      }
      plab <- egf_eval_label(substitute(plab), combined, parent.frame())[subset1]
      if (is.null(plab)) {
        plab <- levels(frame$ts)
      }
      plab <- rep_len(plab, n)
    } else {
      nplot <- n
      main <- egf_eval_label(substitute(main), combined, parent.frame())[subset1]
      if (is.null(main)) {
        main <- levels(frame$ts)
      }
      sub <- egf_eval_label(substitute(sub), combined, parent.frame())[subset1]
      if (is.null(sub)) {
        sub <- ""
      }
      plab <- NULL
    }
    main <- rep_len(main, nplot)
    sub <- rep_len(sub, nplot)

    if (is.null(xlab)) {
      xlab <- switch(time_as, Date = "", numeric = "time")
    }
    if (grepl("^rt(_heat)?$", type)) {
      if (is.null(ylab)) {
        ylab <- "per capita growth rate"
      }
      if (is.null(ylab_outer)) {
        ylab_outer <- "doubling time"
      }
    } else {
      if (is.null(ylab)) {
        ylab <- paste(type, "incidence")
      }
      ylab_outer <- NULL
    }
  }


  ## Augment 'cache' -----------------------------------------------------------

  ## If necessary, initialize
  if (is.null(cache)) {
    cache <- data.frame(
      var = factor(),
      ts = factor(),
      window = factor(),
      time = numeric(0L),
      estimate = numeric(0L),
      se = numeric(0L)
    )
  }
  nc <- names(cache)

  ## If necessary, augment 'cache' with predicted values
  ## of whatever is being plotted and standard errors
  what <- base::sub("_heat$", "", type)
  if (show_predict > 0L) {
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
  if (show_tdoubling > 0L || show_asymptote > 0L) {
    ok <- cache$var == "log(r)" & !(show_tdoubling == 2L & is.na(cache$se))
    required <- setdiff(lw, cache$window[ok])
    if (length(required) > 0L) {
      fx <- fitted(x,
        top = "log(r)",
        link = TRUE,
        se = (show_tdoubling == 2L),
        .subset = (x$frame_windows$window %in% required)
      )
      if (show_tdoubling != 2L) {
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

  ## Augment with confidence intervals
  cache[c("lower", "upper")] <- list(NA_real_)
  i <- (show_predict == 2L & cache$var == sprintf("log(%s)", what)) | (show_tdoubling == 2L & cache$var == "log(r)")
  if (any(i)) {
    cache[i, c("lower", "upper")] <- do.call(do_wald, c(cache[i, c("estimate", "se"), drop = FALSE], list(level = level)))
  }


  ## Augment 'control' ---------------------------------------------------------

  ## Gather parameter values necessary to determine text dimensions ...
  ## those not specified in 'control' are obtained from 'gp'
  par(cex = 1)
  gp <- par()
  if (is.list(control$axis$y)) {
    nel <- setdiff(c("mgp", "cex.axis", "font.axis", "family"), names(control$axis$y))
    control$axis$y[nel] <- gp[nel]
  }
  if (is.list(control$title$sub)) {
    nel <- setdiff(c("mgp", "cex.sub", "font.sub", "family"), names(control$title$sub))
    control$title$sub[nel] <- gp[nel]
  }
  if (is.list(control$title$ylab)) {
    nel <- setdiff(c("mgp", "cex.lab", "font.lab", "family"), names(control$title$ylab))
    control$title$ylab[nel] <- gp[nel]
  }
  if (type == "rt_heat") {
    if (is.list(control$title$plab)) {
      nel <- setdiff(c("cex.lab", "font.lab", "family"), names(control$title$plab))
      control$title$plab[nel] <- gp[nel]
    }
  } else {
    for (s in c("ci", "estimate", "legend")) {
      if (is.list(control$tdoubling[[s]])) {
        nel <- setdiff(c("adj", "cex", "font", "family"), names(control$tdoubling[[s]]))
        control$tdoubling[[s]][nel] <- gp[nel]
      }
    }
  }

  ## Distinguish minor and major axes ...
  if (time_as == "Date" && is.list(args <- control$axis$x)) {
    nel <- setdiff(c("mgp", "cex.axis"), names(args))
    args[nel] <- gp[nel]
    args <- rep_len(list(args), 2L)
    names(args) <- c("minor", "major")
    args$major$mgp[1:2] <- args$minor$mgp[1:2] + 1.5
    args$major$cex.axis <- args$minor$cex.axis * 1.15
    args$major$tick <- FALSE
    control$axis$x <- args
  }


  ## Misc ----------------------------------------------------------------------

  ## Approximation of per capita growth rate from cumulative counts
  ## is much less noisy but sensitive to the (unknown) initial value,
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


  ## Plot and return -----------------------------------------------------------

  if (type == "rt_heat") {
    plot.egf.heat(
      cache = cache,
      time_as = time_as,
      dt = dt,
      log = log,
      per_plot = per_plot,
      control = control,
      xlim = xlim,
      main = main,
      sub = sub,
      xlab = xlab,
      ylab = ylab,
      ylab_outer = ylab_outer,
      plab = plab
    )
  } else {
    plot.egf.curve(
      frame = frame,
      frame_windows = frame_windows,
      cache = cache,
      type = type,
      time_as = time_as,
      dt = dt,
      log = log,
      zero = zero,
      show_predict = show_predict,
      show_tdoubling = show_tdoubling,
      show_asymptote = show_asymptote,
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
plot.egf.curve <- function(frame, frame_windows, cache,
                           type, time_as, dt, log, zero,
                           show_predict, show_tdoubling, show_asymptote,
                           level, control, xlim, ylim,
                           main, sub, xlab, ylab, ylab_outer) {
  n <- nlevels(frame$ts)
  formula <- as.formula(call("~", as.name(type), quote(time)))
  xlim_bak <- xlim
  ylim_bak <- ylim
  xlab_bak <- xlab
  elu <- c("estimate", "lower", "upper")

  ## Graphical parameters
  op <- par(xaxs = "i", yaxs = "i")
  on.exit(par(op))
  gp <- par()

  for (i in seq_len(n)) { # loop over plots
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
        y
      })
    )
    if (log && is.null(zero)) {
      data[[type]][data[[type]] == 0] <- NA
    }
    if (type == "interval") {
      data$pty <- factor(sign(data$dt - dt), levels = c(0, -1, 1), labels = c("main", "short", "long"))
    } else {
      data$pty <- factor("main")
    }

    ## Axis limits (x)
    xlim <- xlim_bak
    if (is.null(xlim)) {
      xlim <- range(data$time)
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
            ylim <- exp(base::log(ry) + c(-1, 1) * 0.04 * diff(base::log(ry)))
            ylim[!is.finite(ylim)] <- ry[!is.finite(ylim)]
          } else {
            if (all(y == 0 | y >= 1) && any(y > 1)) {
              ylim <- exp(c(-0.04, 1.04) * log(max(y)))
            } else {
              ry <- range(y[y > 0])
              ylim <- exp(base::log(ry) + c(-1, 1) * 0.04 * diff(base::log(ry)))
              ylim[!is.finite(ylim)] <- ry[!is.finite(ylim)]
            }
          }
        } else {
          ylim <- c(0, 1.04 * max(y))
        }
      }
    }
    if (log && !is.null(zero)) {
      data[[type]][data[[type]] == 0] <- if (is.na(zero)) ylim[1L] else zero
    }

    plot.new()
    plot.window(xlim = xlim, ylim = ylim, log = if (log) "y" else "")

    ## Fitting windows
    if (is.list(args <- control$window)) {
      args[c("xleft", "xright")] <- data_windows[c("start", "end")]
      args[c("ybottom", "ytop")] <- as.list(ylim)
      do.call(rect, args)
    }

    ## Observed data
    for (s in levels(data$pty)) {
      if (is.list(args <- control$data[[s]])) {
        m <- match(names(args), c("formula", "data", "subset", "x", "y"), 0L)
        args[m] <- NULL
        args <- c(list(formula = formula, data = data, subset = (data$pty == s)), args)
        do.call(points, args)
      }
    }

    ## Confidence bands on predicted curves
    if (show_predict == 2L && is.list(args <- control$predict$ci)) {
      for (px in split(cache_predict, cache_predict$window, drop = TRUE)) {
        if (type == "cumulative") {
          c0 <- data[[type]][match(px$time[1L], data$time, 0L)]
          if (is.na(c0)) {
            next
          }
          px$lower <- c0 + px$lower - px$lower[1L]
          px$upper <- c0 + px$upper - px$upper[1L]
        }
        args$x <- c(px$time, rev(px$time))
        args$y <- c(px$lower, rev(px$upper))
        do.call(polygon, args)
      }
    }

    ## Predicted curves
    if (show_predict > 0L && is.list(args <- control$predict$estimate)) {
      for (px in split(cache_predict, cache_predict$window, drop = TRUE)) {
        if (type == "cumulative") {
          c0 <- data[[type]][match(px$time[1L], data$time, 0L)]
          if (is.na(c0)) {
            next
          }
          px$estimate <- c0 + px$estimate - px$estimate[1L]
        }
        m <- match(names(args), c("formula", "data", "subset", "x", "y"), 0L)
        args[m] <- NULL
        args <- c(list(formula = estimate ~ time, data = px), args)
        do.call(lines, args)
      }
    }

    ## Asymptote
    if (show_asymptote > 0L && is.list(args <- control$asymptote)) {
      args[c("x0", "x1")] <- data_windows[c("start", "end")]
      args[c("y0", "y1")] <- cache_r["estimate"]
      do.call(segments, args)
    }

    ## Box
    if (is.list(args <- control$box)) {
      do.call(box, args)
    }

    ## Axis (x)
    if (is.list(args <- control$axis$x)) {
      args$side <- 1
      do.call(switch(time_as, Date = Daxis, baxis), args)
    }

    ## Axis title (x)
    if (time_as == "numeric" && is.list(args <- control$title$xlab)) {
      args$xlab <- xlab
      do.call(title, args)
    }

    ## Axis (y)
    if (is.list(args1 <- control$axis$y)) {
      args1$side <- 2
      args1$las <- 1
      args1$at <- axTicks(side = 2)
      if (type != "rt" && max(args1$at) >= 1e+05) {
        args1$labels <- get_scientific_labels(args1$at)
      } else {
        args1$labels <- as.character(args1$at)
      }
      do.call(baxis, args1)
      width <- max(strwidth(args1$labels, units = "inches", cex = args1$cex.axis, font = args1$font.axis, family = args1$family))
      line <- args1$mgp[2L] + diff(grconvertX(c(0, width), "inches", "lines")) + 0.75
    }

    ## Axis title (y)
    if (is.list(args2 <- control$title$ylab)) {
      args2$ylab <- ylab
      if (is.list(args1)) {
        args2$line <- line
      } else {
        line <- args2$line
        if (is.null(line)) {
          line <- args2$line <- args2$mgp[1L]
        }
      }
      do.call(title, args2)
      width <- strheight(args2$ylab, units = "inches", cex = args2$cex.lab, font = args2$font.lab, family = args2$family)
      line <- line + diff(grconvertX(c(0, width), "inches", "lines")) + 1.25
    }

    if (type == "rt" && log) {
      ## Axis (y), outer
      if (is.list(args1)) {
        par(new = TRUE)
        plot.window(xlim = xlim, ylim = base::log(2) / ylim, log = "y")
        args1$at <- axTicks(side = 2)
        args1$labels <- as.character(args1$at)
        args1$mgp <- line + args1$mgp
        do.call(baxis, args1)
        width <- max(strwidth(args1$labels, units = "inches", cex = args1$cex.axis, font = args1$font.axis, family = args1$family))
        line <- args1$mgp[2L] + diff(grconvertX(c(0, width), "inches", "lines")) + 0.75
        par(new = TRUE)
        plot.window(xlim = xlim, ylim = ylim, log = "y")
      }

      ## Axis title (y), outer
      if (is.list(args2)) {
        args2$ylab <- ylab_outer
        args2$line <- line
        do.call(title, args2)
      }
    }

    ## Initial doubling times
    if (show_tdoubling > 0L && is.list(args <- control$tdoubling)) {
      tdoubling <- as.list(base::log(2) / cache_r[elu])
      names(tdoubling) <- elu[c(1L, 3L, 2L)]

      tdoubling$labels <- c(
        ci       = sprintf("(%.3g%% CI)", 100 * level),
        estimate = "estimate",
        legend   = "initial doubling time:"
      )
      w <- function(s, l) {
        if (!is.list(l)) {
          return(NA_real_)
        }
        strwidth(s, units = "user", cex = l$cex * gp$cex, font = l$font, family = l$family)
      }
      tdoubling$widths <- mapply(w, s = tdoubling$labels, l = args[names(tdoubling$labels)])
      h <- function(s, l) {
        if (!is.list(l)) {
          return(NA_real_)
        }
        height <- strheight(s, units = "inches", cex = l$cex * gp$cex, font = l$font, family = l$family)
        diff(grconvertY(c(0, height), "inches", "lines"))
      }
      tdoubling$heights <- mapply(h, s = tdoubling$labels, l = args[names(tdoubling$labels)])

      tdoubling$lines <- rep_len(0.25, 5L)
      names(tdoubling$lines) <- paste0(rep.int(c("", "label_"), c(2L, 3L)), names(tdoubling$labels)[c(1:2, 1:3)])

      if (show_tdoubling == 2L && is.list(args$ci)) {
        tdoubling$lines[-1L] <- tdoubling$lines[["ci"]] + tdoubling$heights[["ci"]] + 0.15
      }
      if (is.list(args$estimate)) {
        tdoubling$lines[-(1:2)] <- tdoubling$lines[["estimate"]] + tdoubling$heights[["estimate"]] + 1
      }
      if (show_tdoubling == 2L && is.list(args$ci)) {
        tdoubling$lines[-(1:3)] <- tdoubling$lines[["label_ci"]] + tdoubling$heights[["ci"]] + 0.15
      }
      if (is.list(args$estimate)) {
        tdoubling$lines[-(1:4)] <- tdoubling$lines[["label_estimate"]] + tdoubling$heights[["estimate"]] + 0.25
      }

      ## Legend
      if (show_legend <- is.list(args$legend)) {
        adj <- args$legend$adj
        x_legend <- xlim[1L] + adj * max(0, xlim[2L] - xlim[1L] - tdoubling$widths[["legend"]])
        x_body <- x_legend + tdoubling$widths[["legend"]] - 0.5 * max(tdoubling$widths[c("estimate", "ci")])

        args$legend$text <- tdoubling$labels[["legend"]]
        args$legend$side <- 3
        args$legend$line <- tdoubling$lines[["label_legend"]]
        args$legend$at <- x_legend
        args$legend$adj <- 0
        args$legend$padj <- 0
        args$legend$cex <- args$legend$cex * gp$cex
        do.call(mtext, args$legend)
      }

      ## Estimates
      if (is.list(args$estimate)) {
        args$estimate$text <- c(sprintf("%.1f", tdoubling$estimate), if (show_legend) tdoubling$labels[["estimate"]])
        args$estimate$side <- 3
        args$estimate$line <- c(rep_len(tdoubling$lines[["estimate"]], nrow(data_windows)), if (show_legend) tdoubling$lines[["label_estimate"]])
        args$estimate$at <- c((data_windows$start + data_windows$end) / 2, if (show_legend) x_body)
        args$estimate$adj <- 0.5
        args$estimate$padj <- 0
        args$estimate$cex <- args$estimate$cex * gp$cex
        do.call(mtext, args$estimate)
      }

      ## Confidence intervals
      if (show_tdoubling == 2L && is.list(args$ci)) {
        args$ci$text <- c(sprintf("(%.1f, %.1f)", tdoubling$lower, tdoubling$upper), if (show_legend) tdoubling$labels[["ci"]])
        args$ci$side <- 3
        args$ci$line <- c(rep_len(tdoubling$lines[["ci"]], nrow(data_windows)), if (show_legend) tdoubling$lines[["label_ci"]])
        args$ci$at <- c((data_windows$start + data_windows$end) / 2, if (show_legend) x_body)
        args$ci$adj <- 0.5
        args$ci$padj <- 0
        args$ci$cex <- args$ci$cex * gp$cex
        do.call(mtext, args$ci)
      }
    }

    ## Plot subtitle
    if (is.list(args1 <- control$title$sub)) {
      args1$main <- sub[i]
      names(args1) <- base::sub("\\.sub$", ".main", names(args1))
      if (show_tdoubling > 0L && is.list(control$tdoubling)) {
        args1$line <- tdoubling$lines[["estimate"]] + tdoubling$heights[["estimate"]] + 1
      } else if (is.null(args1$line)) {
        args1$line <- args1$mgp[1L] + 1
      }
      do.call(title, args1)
      height <- strheight(args1$main, units = "inches", cex = args1$cex.main, font = args1$font.main, family = args1$family)
      line <- args1$line + diff(grconvertY(c(0, height), "inches", "lines")) + 0.25
    }

    ## Plot title
    if (is.list(args2 <- control$title$main)) {
      args2$main <- main[i]
      if (is.list(args1)) {
        args2$line <- line
      } else if (show_tdoubling > 0L && is.list(control$tdoubling)) {
        args2$line <- tdoubling$lines[["estimate"]] + tdoubling$heights[["estimate"]] + 1
      }
      do.call(title, args2)
    }
  } # loop over plots

  invisible(NULL)
}

#' @import graphics
#' @importFrom grDevices colorRamp rgb
plot.egf.heat <- function(cache, time_as, dt, log, per_plot, control,
                          xlim, main, sub, xlab, ylab, ylab_outer, plab) {
  n <- nlevels(cache$ts)
  pal <- do.call(colorRamp, control$heat$pal)
  range_log_r <- range(cache$estimate, na.rm = TRUE)
  range_r <- exp(range_log_r)
  if (log) {
    diff_range_log_r <- diff(range_log_r)
    to_unit <- function(x) (x - range_log_r[1L]) / diff_range_log_r
  } else {
    to_unit <- function(x) exp(x) / range_r[2L]
  }
  if (is.null(xlim)) {
    xlim <- range(cache$time)
  }
  xlim1 <- xlim
  ylim1 <- c(0, 1)
  xlim2 <- c(0, 1)
  ylim2 <- if (log) range_r else c(0, range_r[2L])

  ## Graphical parameters
  op <- par(mfrow = c(1, 1), xaxs = "i", yaxs = "i")
  on.exit(par(op))
  mar <- par("mar")

  ## Device layout
  L <- c(seq_len(per_plot), rep_len(per_plot + 1L, per_plot))
  dim(L) <- c(per_plot, 2L)
  layout(L, widths = c(0.925, 0.075))
  par(cex = 1)

  i <- 0L
  while (i < n) { # loop over pages
    for (k in i + seq_len(min(per_plot, n - i))) { # loop over plots
      data <- cache[unclass(cache$ts) == k, , drop = FALSE]

      plot.new()
      plot.window(xlim = xlim1, ylim = ylim1)

      if (k == i + 1L) {
        ## Plot subtitle
        if (is.list(args1 <- control$title$sub)) {
          args1$main <- sub[1 + i / per_plot]
          names(args1) <- base::sub("\\.sub$", ".main", names(args1))
          if (is.null(args1$line)) {
            args1$line <- args1$mgp[1L] + 1
          }
          do.call(title, args1)
          height <- strheight(args1$main, units = "inches", cex = args1$cex.main, font = args1$font.main, family = args1$family)
          line <- args1$line + diff(grconvertY(c(0, height), "inches", "lines")) + 0.25
        }

        ## Plot title
        if (is.list(args2 <- control$title$main)) {
          args2$main <- main[1 + i / per_plot]
          if (is.list(args1)) {
            args2$line <- line
          }
          do.call(title, args2)
        }
      }

      ## Background
      if (is.list(args <- control$heat$bg)) {
        args[c("xleft", "xright", "ybottom", "ytop")] <- as.list(c(xlim1, ylim1))
        do.call(rect, args)
      }

      ## Pixels
      rect(
        xleft = data$time - 0.5 * dt,
        xright = data$time + 0.5 * dt,
        ybottom = ylim1[1L],
        ytop = ylim1[2L],
        border = NA,
        col = rgb(pal(to_unit(data$estimate)), maxColorValue = 255)
      )

      ## Panel title
      if (is.list(args1 <- control$title$plab)) {
        names(args1) <- base::sub("\\.lab$", "", names(args1))

        ## Pad from top left corner
        p <- diff(grconvertY(c(0, 0.08), "npc", "inches"))
        px <- diff(grconvertX(c(0, p), "inches", "user"))
        py <- diff(grconvertY(c(0, p), "inches", "user"))

        ## Underlay
        if (is.list(args2 <- control$heat$ul)) {
          width  <- strwidth( plab[k], units = "user", cex = args1$cex, font = args1$font, family = args1$family)
          height <- strheight(plab[k], units = "user", cex = args1$cex, font = args1$font, family = args1$family)
          args2$xleft   <- xlim1[1L]
          args2$xright  <- xlim1[1L] + width  + 2 * px
          args2$ybottom <- ylim1[2L] - height - 2 * py
          args2$ytop    <- ylim1[2L]
          do.call(rect, args2)
        }

        ## Text
        args1$x <- xlim1[1L] + px
        args1$y <- ylim1[2L] - py
        args1$labels <- plab[k]
        args1$adj <- c(0, 1)
        do.call(text, args1)
      }
    } # loop over plots

    ## Axis (x)
    if (is.list(args <- control$axis$x)) {
      args$side <- 1
      do.call(switch(time_as, Date = Daxis, baxis), args)
    }

    ## Axis title (x)
    if (time_as == "numeric" && is.list(args <- control$title$xlab)) {
      args$xlab <- xlab
      do.call(title, args)
    }

    while (k %% per_plot > 0L) {
      plot.new()
      k <- k + 1L
    }
    par(mar = replace(mar, c(2L, 4L), 0))
    plot.new()
    plot.window(xlim = xlim2, ylim = c(0, 1))

    ## Color scale
    dy <- 0.001
    y <- seq.int(0, 1, by = 2 * dy)
    rect(
      xleft = xlim2[1L],
      xright = xlim2[2L],
      ybottom = y - dy,
      ytop = y + dy,
      border = NA,
      col = rgb(pal(y), maxColorValue = 255)
    )

    ## Axis (y)
    if (is.list(args1 <- control$axis$y)) {
      par(new = TRUE)
      plot.window(xlim = xlim2, ylim = ylim2, log = "y")
      args1$side <- 4
      args1$las <- 1
      args1$at <- axTicks(side = 4)
      args1$labels <- as.character(args1$at)
      do.call(baxis, args1)
      width <- max(strwidth(args1$labels, units = "inches", cex = args1$cex.axis, font = args1$font.axis, family = args1$family))
      line <- args1$mgp[2L] + diff(grconvertX(c(0, width), "inches", "lines")) + 0.75
    }

    ## Axis title (y)
    if (is.list(args2 <- control$title$ylab)) {
      if (!is.list(args1)) {
        line <- args2$line
        if (is.null(line)) {
          line <- args2$mgp[1L]
        }
      }
      tmp <- args2[c("cex.lab", "col.lab", "font.lab", "family", "xpd")]
      names(tmp) <- sub("\\.lab$", "", names(tmp))
      args2 <- list(
        x = xlim2[2L] + diff(grconvertX(c(0, line), "lines", "user")),
        y = if (log) exp(mean(log(ylim2))) else mean(ylim2),
        labels = ylab,
        adj = c(0.5, 0),
        srt = 270
      )
      args2 <- c(args2, tmp)
      do.call(text, args2)
      width <- strheight(ylab, units = "inches", cex = args2$cex, font = args2$font, family = args2$family)
      line <- line + diff(grconvertX(c(0, width), "inches", "lines")) + 1.25
    }

    if (log) {
      ## Axis (y), outer
      if (is.list(args1)) {
        par(new = TRUE)
        plot.window(xlim = xlim2, ylim = base::log(2) / ylim2, log = "y")
        args1$at <- axTicks(side = 4)
        args1$labels <- as.character(args1$at)
        args1$mgp <- line + args1$mgp
        do.call(baxis, args1)
        width <- max(strwidth(args1$labels, units = "inches", cex = args1$cex.axis, font = args1$font.axis, family = args1$family))
        line <- args1$mgp[2L] + diff(grconvertX(c(0, width), "inches", "lines")) + 0.75
        par(new = TRUE)
        plot.window(xlim = xlim2, ylim = ylim2, log = "y")
      }

      ## Axis title (y), outer
      if (is.list(args2)) {
        args2$x = xlim2[2L] + diff(grconvertX(c(0, line), "lines", "user"))
        args2$labels <- ylab_outer
        do.call(text, args2)
      }
    }

    par(mar = mar)
    i <- i + per_plot
  } # loop over pages

  invisible(NULL)
}

