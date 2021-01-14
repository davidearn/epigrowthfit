#' Plot models of epidemic growth
#'
#' @description
#' Method for plotting `"egf"` objects.
#'
#' @param x
#'   An `"egf"` object returned by [egf()].
#' @param subset
#'   An named list of atomic vectors with elements specifying levels
#'   of factors in `object$frame` (and thus fitting windows). Use the
#'   default (`NULL`) if there are no factors in `object$frame` or,
#'   more generally, to plot every fitting window. See Details.
#' @param group_by
#'   A formula of the form `~1`


#' @import graphics
#' @importFrom stats reformulate
plot.egf <- function(x,
                     subset = NULL,
                     group_by = ~1,
                     type = c("interval", "cumulative"),
                     xty = c("date", "numeric"),
                     log = TRUE,
                     legend = FALSE,
                     tol = 0,
                     control = get_control_default("plot.egf"),
                     ...) {

  ### Argument validation #################################

  ## Reduced frame
  frame_red <- x$frame[!duplicated(x$index), -(1:2), drop = FALSE]
  any_factors <- (length(frame_red) > 0L)

  if (any_factors) {
    stop_if_not(
      inherits(group_by, "formula"),
      length(group_by) == 2L,
      grepl("^(1|([[:alnum:]._]+(:[[:alnum:]._]+)*))$", deparse(group_by[[2L]])),
      all.vars(group_by) %in% names(frame_red),
      m = paste0(
        "`group_by` must be a formula of the form `~1` or\n",
        "`~f1:...:fn` with `f1`,...,`fn` factors in `x$frame`."
      )
    )
    group_by <- all.vars(group_by)
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
        m = "`subset` must specify levels of factors in `object$frame`."
      )
      w <- Reduce("&", Map("%in%", frame_red[names(subset)[i]], subset[i]))
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
  xty <- match.arg(xty)
  stop_if_not_tf(log)
  stop_if_not_tf(legend)
  if (type == "interval") {
    stop_if_not(
      is.numeric(tol),
      length(tol) == 1L,
      tol >= 0,
      m = "`tol` must be a non-negative number or `Inf`."
    )
  }

  control_default <- get_control_default("plot.egf")
  if (inherits(control, "list") && !is.null(names(control))) {
    for (pe in intersect(names(control), names(control_default))) {
      if (is.null(control[[pe]])) {
        control_default[pe] <- list(NULL)
      } else if (inherits(control[[pe]], "list") && !is.null(names(control[[pe]]))) {
        s <- intersect(names(control[[pe]]), names(control_default[[pe]]))
        control_default[[pe]][s] <- control[[pe]][s]
      }
    }
  }
  control <- control_default
  dots <- list(...)

  ### Setting up loop over plots ##########################

  ## Augmented frame
  frame_aug <- rbind(x$frame, attr(x$frame, "extra"))
  frame_aug <- cbind(
    index = `length<-`(x$index, nrow(frame_aug)),
    frame_aug
  )

  interaction0 <- function(...) interaction(..., drop = TRUE, sep = ":")
  days <- function(date, since) as.integer(date - since)
  date_diff <- function(date) as.integer(diff(date))

  if (any_factors) {
    ## Split the augmented frame by interaction
    frame_aug_split <- split(frame_aug, interaction0(frame_aug[-(1:3)]))
    ## Split the reduced frame by group of interactions
    if (length(group_by) > 0L) {
      frame_red_split <- split(frame_red, interaction0(frame_red[group_by]))
    } else {
      frame_red_split <- list("1" = frame_red)
    }
    ## Merge grouped interactions
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
      "Plots of multiple time series in one panel are not\n",
      "supported (yet). See `group_by` argument details in\n",
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
  ci <- confint(x, parm = "tdoubling", level = 0.95, method = "wald")

  ### Loop over plots #####################################

  for (d in frame_aug_split) {

    ### Setting up plot ===================================

    index <- droplevels(factor(d[[1L]], levels = index_levels_for_plot))
    date <- d[[2L]]
    cases <- d[[3L]]
    i12 <- vapply(levels(index), function(s) range(which(index == s])), integer(2L))
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
      dt = c(NA, date_diff(date))
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
      pr <- predict(x,
        subset = if (length(d) > 3L) d[-(1:3), i1, drop = FALSE],
        time = seq.int(from = 0L, to = t2 - t1, by = by),
        se = TRUE
      )
      ci <- confint(pred, level = 0.95, log = FALSE)[[varname]]
      if (type == "cumulative" && i1 > 1L) {
        ci[, -1L] <- sum(cases[2L:i1]) + ci[, -1L]
      }
      ci[[1L]] <- t1 + ci[[1L]]
      ci
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
      if (any_groups) {
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

    op <- par(
      mar = c(4, 5, 2.7, 0.5 + 6 * annotate) + 0.1,
      bty = "l",
      xaxs = "i",
      yaxs = "i",
      las = 1
    )
    on.exit({
      assign("egf.par", par("mar", "plt"), envir = .epigrowthfit)
      par(op)
    })

    plot.new()
    plot.window(xlim = xlim, ylim = ylim, log = if (log) "y" else "")

    ## Fitting windows
    if (!is.null(control$window)) {
      for (i in seq_len(nlevels(index))) {
        l <- list(
          x = c(t1[i], t2[i], t2[i], t1[i]),
          y = ylim[c(1, 1, 2, 2)]
        )
        do.call(polygon, c(l, control$window))
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
        l <- list(side = 1L)
        control$xax$mgp <- c(3, control$xax$mgp2, 0)
        control$xax$mgp2 <- NULL
        do.call(axis, c(l, control$xax))
      }
    }

    ## Axis (y)
    if (!is.null(control$yax)) {
      yax_at <- axTicks(side = 2L)
      if (max(yax_at) < 1e05) {
        yax_labels <- TRUE
        long_yax_labels_flag <- FALSE
      } else {
        yax_labels <- get_labels(yax_at)
        mlw <- max(strwidth(yax_labels, units = "inches", cex = 0.85))
        long_yax_labels_flag <- (mlw / par("csi") + control$yax$mgp2 > 3.75)
      }
      l <- list(
        side = 2L,
        at = yax_at,
        labels = yax_labels
      )
      control$yax$mgp <- c(3, control$yax$mgp2, 0)
      control$yax$mgp2 <- NULL
      if (is.na(control$yax$cex.axis)) {
        control$yax$cex.axis <- (1 - 0.25 * long_yax_labels_flag) * 0.85
      }
      do.call(axis, c(l, control$yax))
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
    if (!is.null(control$confband)) {
      for (p in pred) {
        l <- list(
          x = c(p$time, rev(p$time)),
          y = c(p$lower, rev(p$upper))
        )
        do.call(polygon, c(l, control$confband))
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
    if (!is.null(control$text_dbl) && type == "interval" && x$curve %in% c("exponential", "logistic", "richards")) {
      for (s in levels(index)) {
        elu <- unlist(ci[match(s, levels(x$index)), length(ci) - 2:0])
        l <- list(
          labels = sprintf("doubling time:\n%.1f (%.1f, %.1f) days", elu[1L], elu[2L], elu[3L]),
          xpd = NA
        )
        if (is.na(control$text_dbl$y)) {
          er <- range(pred[[s]]$estimate, na.rm = TRUE)
          if (log) {
            control$text_dbl$y <- er[1L] * 10^(0.25 * diff(log10(er)))
          } else {
            control$text_dbl$y <- er[1] + 0.25 * diff(er)
          }
        }
        if (is.na(control$text_dbl$x)) {
          control$text_dbl$x <- min(pred[[s]]$time[pred[[s]]$estimate > control$text_dbl$y], na.rm = TRUE)
        } else {
          if (is.character(control$text_dbl$x)) {
            control$text_dbl$x <- as.Date(control$text_dbl$x)
          }
          if (inherits(control$text_dbl$x, "Date")) {
            control$text_dbl$x <- days(control$text_dbl$x, since = d0)
          }
        }
        do.call(text, c(l, control$text_dbl))
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
        pch = c(get_pe("pch"), NA)[lind],
        pt.bg = c(get_pe("bg"), NA)[index],
        lty = c(rep(NA, 3L), control$lines$lty)[lind],
        lwd = c(rep(NA, 3L), control$lines$lwd)[lind],
        col = c(get_pe("col"), control$lines$col)[lind]
      )
    }
  }

  invisible(NULL)
}
