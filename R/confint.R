#' Confidence intervals on fitted values
#'
#' Computes confidence intervals on fitted values (see [fitted.egf()])
#' of nonlinear model parameters and, where appropriate,
#' the basic reproduction number and initial doubling time (in days).
#'
#' @param object
#'   An `"egf"` object returned by [egf()].
#' @param parm
#'   Unused argument included for generic consistency.
#' @param level
#'   A number in the interval (0,1). The desired confidence level.
#' @param par
#'   A subset of `get_par_names(object, link = TRUE)` naming nonlinear
#'   model parameters for which confidence intervals should be computed.
#'   If `object$curve %in% c("exponential", "logistic", "richards")`,
#'   then `par` may also contain `"R0"` and `"tdoubling"`.
#' @param subset
#'   An expression to be evaluated in the combined model frame
#'   (see [make_combined()]). Must evaluate to a logical vector
#'   indexing rows of the data frame, and thus fitting windows.
#'   Confidence intervals are computed only for indexed windows.
#'   The default (`NULL`) is to consider all windows.
#' @param link
#'   A logical scalar. If `FALSE`, then confidence intervals
#'   on inverse link-transformed fitted values are returned.
#' @param method
#'   A character string indicating how confidence intervals
#'   should be calculated (see Details).
#' @param grid_len (For `method = "profile"`.)
#'   A positive integer. Step sizes chosen adaptively by
#'   [TMB::tmbprofile()] will generate approximately this
#'   many points on each side of a profile's minimum point.
#' @param max_width (For `method = "uniroot"`.)
#'   A positive number. [TMB::tmbroot()] will search for roots
#'   in the interval from `x-max_width` to `x+max_width`, where
#'   `x` is the fitted value (link scale).
#' @param trace (For `method != "wald"`.)
#'   A logical scalar. If `TRUE`, then basic tracing messages
#'   are printed.
#' @param breaks,probs (For `parm = "R0"`.)
#'   Arguments to `compute_R0()`.
#' @inheritParams check_parallel
#' @inheritParams fitted.egf
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' Three methods are provided for calculating confidence intervals:
#' \describe{
#' \item{`wald`}{
#'   See [fitted.egf()] and [confint.egf_fitted()].
#' }
#' \item{`profile`}{
#'   See [profile.egf()] and [confint.egf_profile()].
#' }
#' \item{`uniroot`}{
#'   Similar to `"profile"`, except that the two solutions of
#'   `deviance(value) = qchisq(level, df = 1)` are approximated
#'   by root-finding using [TMB::tmbroot()] ([stats::uniroot()]
#'   internally).
#' }
#' }
#' For nonlinear model parameters following random effects models,
#' `"wald"` returns confidence intervals on individual fitted values,
#' whereas `"profile"` and `"uniroot"` return confidence intervals
#' on the fixed effects components only, namely the population fitted
#' values.
#'
#' `"wald"` requires minimal computation time but assumes, e.g.,
#' asymptotic normality of the maximum likelihood estimator.
#' A further limitation of `"wald"` is functional non-invariance.
#' `"profile"` and `"uniroot"` avoid these issues but are much
#' slower, requiring estimation of restricted models. Of the two,
#' `"profile"` is more robust.
#'
#' See topic [`nse`] for details on nonstandard evaluation
#' of `subset` and `append`.
#'
#' @return
#' A data frame inheriting from class `"egf_confint"`,
#' with variables:
#' \item{`par`}{
#'   Nonlinear model parameter,
#'   from `get_par_names(object, link = TRUE)`.
#' }
#' \item{`ts`}{
#'   Time series, from `levels(object$endpoints$ts)`.
#' }
#' \item{`window`}{
#'   Fitting window, from `levels(object$endpoints$window)`.
#' }
#' \item{`estimate`, `lower`, `upper`}{
#'   Fitted value and approximate lower and upper confidence limits.
#' }
#' `level` and `object$endpoints` are retained as attributes.
#'
#' @seealso [plot.egf_confint()]
#' @export
#' @importFrom stats qchisq fitted profile confint
#' @importFrom TMB tmbroot
#' @import parallel
confint.egf <- function(object,
                        parm,
                        level = 0.95,
                        par = get_par_names(object, link = TRUE),
                        subset = NULL,
                        append = NULL,
                        link = TRUE,
                        method = c("wald", "profile", "uniroot"),
                        grid_len = 12,
                        max_width = 7,
                        trace = TRUE,
                        breaks = NULL,
                        probs = NULL,
                        parallel = c("serial", "multicore", "snow"),
                        cores = getOption("egf.cores", 2L),
                        outfile = NULL,
                        cl = NULL,
                        .subset = NULL,
                        .append = NULL,
                        ...) {
  stop_if_not_number_in_interval(level, 0, 1, "()")
  stop_if_not_true_false(link)
  method <- match.arg(method)
  s_elu <- c("estimate", "lower", "upper")

  pn <- pn_bak <- get_par_names(object, link = TRUE)
  spec <- c("R0", "tdoubling")
  if (object$curve %in% c("exponential", "logistic", "richards")) {
    pn <- c(pn, spec)
  }

  par <- par_bak <- unique(match.arg(par, pn, several.ok = TRUE))
  if ("R0" %in% par) {
    stop_if_not(
      !is.null(breaks),
      !is.null(probs),
      m = "par = \"R0\": `breaks` and `probs` must be non-NULL.\nSee `help(\"compute_R0\", \"epigrowthfit\")`."
    )
  }
  par[par %in% spec] <- "log(r)"
  par <- unique(par)

  combined <- make_combined(object)
  subset <- subset_to_index(substitute(subset), combined, parent.frame(),
                            .subset = .subset)
  .subset <- replace(rep_len(FALSE, nrow(combined)), subset, TRUE)
  append <- append_to_index(substitute(append), combined, parent.frame(),
                            .append = .append)
  .append <- append

  if (method == "wald") {
    ft <- fitted(object, par = par, se = TRUE, .subset = .subset, .append = .append)
    out <- confint(ft, level = level, link = link)

  } else if (method == "profile") {
    pf <- profile(object, par = par, .subset = .subset, .append = .append,
      max_level = level + min(0.01, 0.1 * (1 - level)),
      grid_len = grid_len,
      trace = trace,
      parallel = parallel,
      cores = cores,
      outfile = outfile,
      cl = cl
    )
    out <- confint(pf, level = level, link = link)
    out$linear_combination <- NULL

  } else { "uniroot"
    stop_if_not_number_in_interval(max_width, 0, Inf, "()")

    p <- length(par)
    w <- length(subset)
    m <- p * w
    n <- length(object$nonrandom)
    A <- rep.int(0, m * n)
    dim(A) <- c(m, n)
    for (k in seq_len(p)) {
      i <- (k - 1L) * w + seq_len(w)
      j <- which(object$tmb_args$data$X_info$par == par[k])
      A[i, j] <- object$tmb_args$data$X[subset, j]
    }
    A_rows <- lapply(seq_len(m), function(i) A[i, ])

    a <- list(
      obj = object$tmb_out,
      target = qchisq(level, df = 1) / 2, # y := diff(nll) = deviance / 2
      sd.range = max_width,
      trace = FALSE
    )
    i_of_m <- function(i) {
      sprintf("%*d of %d", nchar(m), i, m)
    }

    do_uniroot <- function(r, i) {
      if (trace) {
        cat("Computing confidence interval", i_of_m(i), "...\n")
      }
      a$lincomb <- r
      do.call(tmbroot, a)
    }

    if (parallel == "snow") {
      environment(do_uniroot) <- environment(i_of_m) <- .GlobalEnv # see R/boot.R
      if (is.null(cl)) {
        if (is.null(outfile)) {
          outfile <- ""
        }
        cl <- makePSOCKcluster(cores, outfile = outfile)
        on.exit(stopCluster(cl))
      }
      clusterEvalQ(cl, library("TMB"))
      clusterExport(cl, varlist = c("a", "i_of_m", "m"), envir = environment())
      lul <- clusterMap(cl, do_uniroot, r = A_rows, i = seq_len(m))
    } else {
      if (!is.null(outfile)) {
        sink(outfile, type = "output")
        sink(outfile, type = "message")
      }
      if (parallel == "multicore") {
        lul <- mcmapply(do_uniroot, r = A_rows, i = seq_len(m),
                        SIMPLIFY = FALSE, mc.cores = cores)
      } else { # serial
        lul <- Map(do_uniroot, r = A_rows, i = seq_len(m))
      }
      if (!is.null(outfile)) {
        sink(type = "output")
        sink(type = "message")
      }
    }

    out <- data.frame(
      par = rep(factor(par, levels = pn), each = w),
      ts = object$endpoints$ts[subset],
      window = object$endpoints$window[subset],
      estimate = as.numeric(A %*% object$best[object$nonrandom]),
      do.call(rbind, lul), # "lower" "upper"
      combined[subset, append, drop = FALSE],
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    if (!link) {
      out[s_elu] <- apply_inverse_link(out[s_elu], g = out$par)
      levels(out$par) <- remove_link_string(levels(out$par))
    }
  }

  if (any(par_bak %in% spec)) {
    s <- if (link) "log(r)" else "r"
    d_r <- out[out$par == s, , drop = FALSE]
    if (link) {
      d_r[s_elu] <- exp(d_r[s_elu])
    }

    if ("R0" %in% par_bak) {
      d_R0 <- d_r
      d_R0[s_elu] <- lapply(d_r[s_elu], compute_R0, breaks = breaks, probs = probs)
      d_R0$par <- factor("R0")
      out <- rbind(out, d_R0)
    }
    if ("tdoubling" %in% par_bak) {
      s_eul <- s_elu[c(1L, 3L, 2L)]
      d_tdoubling <- d_r
      d_tdoubling[s_elu] <- log(2) / d_r[s_eul]
      d_tdoubling$par <- factor("tdoubling")
      out <- rbind(out, d_tdoubling)
    }
    if (!"log_r" %in% par_bak) {
      out <- out[out$par != s, , drop = FALSE]
    }
  }

  row.names(out) <- NULL
  attr(out, "level") <- level
  attr(out, "endpoints") <- object$endpoints # for `plot.egf_confint()`
  class(out) <- c("egf_confint", "data.frame")
  out
}

#' Plot confidence intervals
#'
#' A method for graphically comparing confidence intervals
#' on fitted values of nonlinear model parameters across
#' fitting windows.
#'
#' @param x
#'   An `"egf_confint"` object returned by [confint.egf()].
#' @param type
#'   A character string determining how confidence intervals
#'   are displayed (see Details).
#' @param subset
#'   An expression to be evaluated in `x`. Must evaluate to
#'   a logical vector indexing rows of `x`. Only indexed
#'   confidence intervals are plotted. The default (`NULL`)
#'   is to plot all confidence intervals.
#' @param order
#'   An expression to be evaluated in `x`, typically a call
#'   to [order()], determining the order in which confidence
#'   intervals or time series (depending on `type`) are plotted.
#'   Must evaluate to a permutation of `seq_len(nrow(x))`.
#'   The default (`NULL`) is equivalent to `seq_len(nrow(x))`.
#' @param per_plot
#'   A positive integer. One plot will display at most this many
#'   confidence intervals or time series (depending on `type`).
#' @param main
#'   An expression or character string indicating a plot title,
#'   to be recycled for all plots.
#' @param label
#'   An expression to be evaluated in `x`, typically a factor,
#'   interaction of factor, or call to [sprintf()]. Used to
#'   create appropriate _y_-axis labels for confidence intervals
#'   when `type = "bars"` and approprate panel labels for
#'   time series when `type = "boxes"`. The default (`NULL`)
#'   is to take labels from `x$window` and `x$ts`, respectively.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' `type = "bars"` creates a one-dimensional plot with
#' confidence intervals drawn as stacked horizontal line segments.
#'
#' `type = "boxes"` creates a two-dimensional plot for each
#' time series (level of `x$ts`), each containing shaded boxes.
#' Projection of boxes onto the horizontal and vertical axes
#' yields fitting windows and corresponding confidence intervals,
#' respectively.
#'
#' If an endpoint of a confidence interval is `NA`, then dashed
#' lines are drawn from the point estimate to the boundary of the
#' plotting region to indicate missingness.
#'
#' @return
#' `NULL` (invisibly).
#'
#' @export
plot.egf_confint <- function(x,
                             type = c("bars", "boxes"),
                             subset = NULL,
                             order = NULL,
                             per_plot = switch(type, bars = 12L, 4L),
                             main = NULL,
                             label = NULL,
                             ...) {
  type <- match.arg(type)
  stop_if_not_positive_integer(per_plot)

  subset <- subset_to_index(substitute(subset), x, parent.frame())
  order <- order_to_index(substitute(order), x, parent.frame())

  a <- attributes(x)
  x <- x[order[order %in% subset], , drop = FALSE]

  if (is.null(main)) {
    s <- switch(type, bars = "fitting window", "time series")
    main <- sprintf("%.3g%% CI by %s", 100 * a$level, s)
  }

  label <- label_to_character(substitute(label), x, parent.frame())
  if (is.null(label)) {
    s <- switch(type, bars = "window", "ts")
    label <- as.character(x[[s]])
  }

  nx <- c("par", "ts", "window", "estimate", "lower", "upper")
  x <- data.frame(x[nx], label, stringsAsFactors = FALSE)

  if (type == "bars") {
    do_bars_plot(x, per_plot = per_plot, main = main)
  } else {
    shift <- min(a$endpoints$start)
    origin <- attr(a$endpoints, "origin") + shift
    iep <- match(x$window, a$endpoints$window, 0L)
    endpoints <- a$endpoints[iep, c("start", "end"), drop = FALSE] - shift
    x <- data.frame(x, endpoints)

    do_boxes_plot(x, origin = origin, per_plot = per_plot, main = main)
  }
  invisible(NULL)
}

#' @keywords internal
#' @import graphics
do_bars_plot <- function(x, per_plot, main) {
  mar <- c(3.5, 5, 1.5, 1)
  csi <- par("csi")
  op <- par(mar = mar)
  on.exit(par(op))

  x_split <- split(x, factor(x$par, levels = unique(x$par)))
  nxs <- names(x_split)

  yax_cex <- get_yax_cex(x$label, mex = 0.92 * mar[2L], cex = 0.8, font = 1, csi = csi)
  yax_cex <- min(0.8, yax_cex)

  for (i in seq_along(x_split)) { # loop over nonlinear model parameters
    xi <- x_split[[i]]
    xlab <- nxs[i]
    xlim <- range(xi[c("estimate", "lower", "upper")], na.rm = TRUE)
    lna <- is.na(xi$lower)
    una <- is.na(xi$upper)

    K <- 0L
    while (K < nrow(xi)) { # loop over plots
      k <- K + seq_len(min(per_plot, nrow(xi) - K))
      plot.new()
      plot.window(xlim = xlim, ylim = c(per_plot + 1, 0), xaxs = "r", yaxs = "i")
      usr <- par("usr")
      abline(v = axTicks(side = 1), lty = 3, col = "grey75")
      segments(
        x0  = replace(xi$lower[k], lna[k], usr[1L]),
        x1  = xi$estimate[k],
        y0  = seq_along(k),
        y1  = seq_along(k),
        lty = c(1, 2)[1L + lna[k]],
        lwd = c(2, 1)[1L + lna[k]]
      )
      segments(
        x0  = xi$estimate[k],
        x1  = replace(xi$upper[k], una[k], usr[2L]),
        y0  = seq_along(k),
        y1  = seq_along(k),
        lty = c(1, 2)[1L + una[k]],
        lwd = c(2, 1)[1L + una[k]]
      )
      points(
        x   = xi$estimate[k],
        y   = seq_along(k),
        pch = 21,
        bg  = "grey75"
      )
      box()
      axis(
        side = 1,
        tcl = -0.4,
        mgp = c(3, 0.3, 0),
        cex.axis = yax_cex
      )
      axis(
        side = 2,
        at = seq_along(k),
        labels = xi$label[k],
        tick = FALSE,
        las = 1,
        mgp = c(3, 0.2, 0),
        cex.axis = 0.8
      )
      title(xlab = xlab, line = 2, cex.lab = 0.9)
      title(main, line = 0.25, adj = 0, cex.main = 0.9)
      K <- K + per_plot
    } # loop over plots
  } # loop over nonlinear model parameters

  invisible(NULL)
}

#' @keywords internal
#' @import graphics
do_boxes_plot <- function(x, origin, per_plot, main) {
  op <- par(
    mfrow = c(per_plot, 1),
    mar = c(0, 3, 0.25, 0.5),
    oma = c(2.5, 1.5, 1.5, 0)
  )
  on.exit(par(op))

  x_split <- split(x, factor(x$par, levels = unique(x$par)))
  nxs <- names(x_split)

  xlim <- c(0, max(x$end))
  xax_at <- daxis(origin = origin + 1, plot = FALSE)

  for (i in seq_along(x_split)) { # loop over nonlinear model parameters
    xi <- x_split[[i]]
    xi_split <- split(xi, factor(xi$ts, levels = unique(xi$ts)))
    ylab <- nxs[i]
    ylim <- range(xi[c("estimate", "lower", "upper")], na.rm = TRUE)

    K <- 0L
    while (K < length(xi_split)) { # loop over plots
      for (k in K + seq_len(min(per_plot, length(xi_split) - K))) { # loop over panels
        xik <- xi_split[[k]]

        plot.new()
        plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "r")
        gp <- par(c("usr", "cex", "mai", "omi", "pin", "din"))
        abline(v = xax_at, lty = 3, col = "grey75")

        xik$lower[lna <- is.na(xik$lower)] <- gp$usr[3L]
        xik$upper[una <- is.na(xik$upper)] <- gp$usr[4L]

        for (l in seq_len(nrow(xik))) { # loop over boxes
          rect(
            xleft   = xik$start[l],
            xright  = xik$end[l],
            ybottom = xik$lower[l],
            ytop    = xik$estimate[l],
            col     = if (lna[l]) NA else "grey75",
            border  = "grey50",
            lty     = if (lna[l]) 2 else 1
          )
          rect(
            xleft   = xik$start[l],
            xright  = xik$end[l],
            ybottom = xik$estimate[l],
            ytop    = xik$upper[l],
            col     = if (una[l]) NA else "grey75",
            border  = "grey50",
            lty     = if (una[l]) 2 else 1
          )
          segments(
            x0  = xik$start[l],
            x1  = xik$end[l],
            y0  = xik$estimate[l],
            y1  = xik$estimate[l],
            col = "grey50",
            lwd = 2
          )
        } # loop over boxes

        box()
        axis(
          side = 2,
          las = 1,
          mgp = c(3, 0.7, 0),
          cex.axis = 0.8
        )
        text(
          x = gp$usr[1L] + 0.075 * (gp$usr[2L] - gp$usr[1L]) * (gp$pin[2L] / gp$pin[1L]),
          y = gp$usr[4L] - 0.075 * (gp$usr[4L] - gp$usr[3L]),
          labels = xik$label[1L],
          adj = c(0, 1),
          cex = 0.8
        )
      } # loop over panels

      mtext(main,
        side = 3,
        line = 0,
        outer = TRUE,
        at = gp$mai[2L] / (gp$din[1L] - sum(gp$omi[c(2L, 4L)])),
        adj = 0,
        cex = 0.9,
        font = 2
      )
      mtext(ylab,
        side = 2,
        line = 0,
        outer = TRUE,
        at = 0.5,
        adj = 0.5,
        cex = 0.9
      )
      daxis(
        origin = origin + 1,
        minor = list(
          mgp = c(3, 0.25, 0), xpd = TRUE,
          lwd = 0, lwd.ticks = 1, tcl = -0.2,
          gap.axis = 0, cex.axis = 0.8 / gp$cex
        ),
        major = list(
          mgp = c(3, 1.25, 0), xpd = TRUE,
          lwd = 0, lwd.ticks = 0, tcl = 0,
          gap.axis = 0, cex.axis = 0.9 / gp$cex
        )
      )
      K <- K + per_plot
    } # loop over plots
  } # loop over nonlinear model parameters

  invisible(NULL)
}
