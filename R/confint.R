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
#'   An expression to be evaluated in the combined model frame.
#'   Must evaluate to a logical vector or list of logical vectors
#'   indexing rows of the model frame, and thus fitting windows.
#'   Confidence intervals are computed only for the indexed
#'   fitting windows. The default (`NULL`) is to consider all
#'   fitting windows.
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
#' @inheritSection fitted.egf Nonstandard evaluation
#' @inheritSection fitted.egf Warning
#'
#' @return
#' A data frame inheriting from class `"egf_confint"`,
#' with variables:
#' \item{`par`}{
#'   Nonlinear model parameter,
#'   from `get_par_names(object, link = TRUE)`.
#' }
#' \item{`ts`}{
#'   Time series, from `levels(object$endpoints$window)`.
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
  par[par %in% spec] <- "log_r"
  par <- unique(par)

  frame <- do.call(cbind, unname(object$frame_par))
  frame <- frame[unique(names(frame))]
  subset <- subset_to_index(substitute(subset), frame, parent.frame(),
                            .subset = .subset)
  .subset <- replace(rep_len(FALSE, nrow(frame)), subset, TRUE)
  append <- append_to_index(substitute(append), frame, parent.frame(),
                            .append = .append)
  .append <- names(frame[append])

  if (method == "wald") {
    ft <- fitted(object, par = par, .subset = .subset, .append = .append)
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
      frame[subset, append, drop = FALSE],
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
    s <- if (link) "log_r" else "r"
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
  structure(out,
    level = level,
    endpoints = object$endpoints, # for `plot.egf_confint()`
    class = c("egf_confint", "data.frame")
  )
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
#'   An expression to be evaluated in `x`. Must evaluate to a
#'   logical vector or list of logical vectors indexing rows
#'   of `x`. Only indexed confidence intervals are plotted.
#'   The default (`NULL`) is to plot all confidence intervals.
#' @param order
#'   An expression to be evaluated in `x`, typically a call
#'   to [order()], determining the order in which confidence
#'   intervals and time series (depending on `type`) are
#'   plotted. Must evaluate to a permutation of `seq_len(nrow(x))`.
#'   The default (`NULL`) is equivalent to `seq_len(nrow(x))`.
#' @param label
#'   An expression to be evaluated in `x`, typically a factor,
#'   interaction of factors, or call to [sprintf()]. Used to
#'   create appropriate labels for confidence intervals and
#'   time series (depending on `type`). The default (`NULL`)
#'   is to take labels from `x$window` and `x$ts`, respectively.
#' @param per_plot
#'   A positive integer. One plot will display at most this many
#'   confidence intervals or time series (depending on `type`).
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
#' yields fitting windows and corresponding confidence intervals.
#'
#' If an endpoint of a confidence interval is `NA`, then dashed
#' lines are drawn from the point estimate to the boundary of the
#' plotting region to indicate missingness.
#'
#' @return
#' `NULL` (invisibly).
#'
#' @export
#' @import graphics
plot.egf_confint <- function(x,
                             type = c("bars", "boxes"),
                             subset = NULL,
                             order = NULL,
                             label = NULL,
                             per_plot = switch(type, bars = 12L, 4L),
                             ...) {
  type <- match.arg(type)
  stop_if_not_positive_integer(per_plot)

  subset <- subset_to_index(substitute(subset), x, parent.frame())
  order <- order_to_index(substitute(order), x, parent.frame())
  label <- substitute(label)
  if (is.null(label)) {
    label <- as.character(x[[switch(type, bars = "window", "ts")]])
    label_format <- switch(type, bars = "fitting window", "time series")
  } else {
    label <- label_to_character(label, x, parent.frame())
    label_format <- attr(label, "format")
  }
  a <- attributes(x)
  x <- cbind(label, x)
  x <- x[order[subset], , drop = FALSE]
  if (type == "boxes") {
    shift <- min(a$endpoints$start)
    origin <- attr(a$endpoints, "origin") + shift
    x <- cbind(a$endpoints[match(x$window, a$endpoints$window, 0L), c("start", "end"), drop = FALSE] - shift, x)
  }
  x_split <- split(x, x$par, drop = TRUE)

  if (type == "bars") {
    op <- par(mar = c(3.5, 5, 1.5, 1))
    on.exit(par(op))

    msw <- max(strwidth(x$label, units = "inches", cex = 0.8, font = 1))
    yax_cex <- 0.8 / max(1, msw / (0.92 * par("mai")[2L]))

    for (i in seq_along(x_split)) { # loop over nonlinear model parameters
      xi <- x_split[[i]]
      xlab <- names(x_split)[i]
      xlim <- range(xi[c("estimate", "lower", "upper")], na.rm = TRUE)
      is_na_lower <- is.na(xi$lower)
      is_na_upper <- is.na(xi$upper)

      j <- 0L
      while (j < nrow(xi)) { # loop over plots
        k <- j + seq_len(min(per_plot, nrow(xi) - j))
        plot.new()
        plot.window(xlim = xlim, ylim = c(per_plot + 1, 0), yaxs = "i")
        abline(v = axTicks(side = 1), lty = 3, col = "grey75")
        segments(
          x0 = ifelse(is_na_lower[k], par("usr")[1L], xi$lower[k]),
          x1 = xi$estimate[k],
          y0 = seq_along(k),
          y1 = seq_along(k),
          lty = c(1, 2)[1L + is_na_lower[k]],
          lwd = c(2, 1)[1L + is_na_lower[k]]
        )
        segments(
          x0 = xi$estimate[k],
          x1 = ifelse(is_na_upper[k], par("usr")[2L], xi$upper[k]),
          y0 = seq_along(k),
          y1 = seq_along(k),
          lty = c(1, 2)[1L + is_na_upper[k]],
          lwd = c(2, 1)[1L + is_na_upper[k]]
        )
        points(
          x = xi$estimate[k],
          y = seq_along(k),
          pch = 21,
          bg = "grey75"
        )
        box()
        axis(side = 1,
          tcl = -0.4,
          mgp = c(3, 0.3, 0),
          cex.axis = yax_cex
        )
        axis(side = 2,
          at = seq_along(k),
          labels = xi$label[k],
          tick = FALSE,
          las = 1,
          mgp = c(3, 0.2, 0),
          cex.axis = 0.8
        )
        title(xlab = xlab, line = 2, cex.lab = 0.9)
        main <- sprintf("%.3g%% CI by %s", 100 * a$level, label_format)
        title(main, line = 0.25, adj = 0, cex.main = 0.9)
        j <- j + per_plot
      } # loop over plots
    } # loop over nonlinear model parameters

  } else { # "boxes"
    op <- par(
      mfrow = c(per_plot, 1),
      mar = c(0, 3, 0.25, 0.5),
      oma = c(2.5, 1.5, 1.5, 0)
    )
    on.exit(par(op))

    xlim <- c(0L, max(x$end))
    xax_at <- daxis(
      left = xlim[1L],
      right = xlim[2L],
      origin = origin + 1,
      plot = FALSE
    )

    for (i in seq_along(x_split)) { # loop over nonlinear model parameters
      xi <- x_split[[i]]
      ylim <- range(xi[c("estimate", "lower", "upper")], na.rm = TRUE)
      xi_split <- split(xi, xi$ts, drop = TRUE)

      j <- 0L
      while (j < length(xi_split)) { # loop over plots
        for (k in j + seq_len(min(per_plot, length(xi_split) - j))) { # loop over panels
          xik <- xi_split[[k]]
          plot.new()
          plot.window(xlim = xlim, ylim = ylim, xaxs = "i")
          abline(v = xax_at, lty = 3, col = "grey75")

          for (l in seq_len(nrow(xik))) { # loop over boxes
            t12 <- unlist(xik[l, c("start", "end")], use.names = FALSE)
            elu <- unlist(xik[l, c("estimate", "lower", "upper")], use.names = FALSE)
            if (lna <- is.na(elu[2L])) {
              elu[2L] <- par("usr")[3L]
            }
            if (una <- is.na(elu[3L])) {
              elu[3L] <- par("usr")[4L]
            }
            rect(
              xleft = t12[1L],
              xright = t12[2L],
              ybottom = elu[2L],
              ytop = elu[1L],
              col = if (lna) NA else "grey75",
              border = "grey50",
              lty = if (lna) 2 else 1
            )
            rect(
              xleft = t12[1L],
              xright = t12[2L],
              ybottom = elu[1L],
              ytop = elu[3L],
              col = if (una) NA else "grey75",
              border = "grey50",
              lty = if (una) 2 else 1
            )
            segments(
              x0 = t12[1L],
              x1 = t12[2L],
              y0 = elu[1L],
              y1 = elu[1L],
              col = "grey50",
              lwd = 2,
            )
          } # loop over boxes

          box()
          axis(side = 2, las = 1, mgp = c(3, 0.7, 0), cex.axis = 0.85)
          text(
            x = par("usr")[1L] + (0.075 * par("pin")[2L] / par("pin")[1L]) * diff(par("usr")[1:2]),
            y = par("usr")[4L] - 0.075 * diff(par("usr")[3:4]),
            labels = xik$label[1L],
            adj = c(0, 1),
            cex = 0.8
          )
        } # loop over panels

        main <- sprintf("%.3g%% CI by %s", 100 * a$level, label_format)
        mtext(main,
          side = 3,
          line = 0,
          outer = TRUE,
          at = par("mai")[2L] / (par("din")[1L] - sum(par("omi")[c(2L, 4L)])),
          adj = 0,
          cex = 0.8,
          font = 2
        )
        mtext(names(x_split)[i],
          side = 2,
          line = 0,
          outer = TRUE,
          at = 0.5,
          adj = 0.5,
          cex = 0.8
        )
        daxis(
          left = xlim[1L],
          right = xlim[2L],
          origin = origin + 1,
          minor = list(mgp = c(3, 0.25, 0), tcl = -0.2, gap.axis = 0,
                       lwd = 0, lwd.ticks = 1, cex.axis = 0.85, xpd = TRUE),
          major = list(mgp = c(3, 1.25, 0), tcl = 0,    gap.axis = 0,
                       lwd = 0, lwd.ticks = 0, cex.axis = 0.85, xpd = TRUE)
        )
        j <- j + per_plot
      } # loop over plots
    } # loop over nonlinear model parameters
  }

  invisible(NULL)
}
