#' Confidence intervals on fitted values
#'
#' Computes confidence intervals on \link[=fitted.egf]{fitted values}
#' of top level nonlinear model parameters and, where appropriate,
#' initial doubling times and basic reproduction numbers.
#'
#' @param object
#'   An \code{"\link{egf}"} object.
#' @param parm
#'   Unused argument included for generic consistency.
#' @param level
#'   A number in the interval (0,1) indicating a confidence level.
#' @param top
#'   A subset of \code{\link{egf_get_names_top}(object, link = TRUE)}
#'   naming top level nonlinear model parameters for which confidence
#'   intervals should be computed. If \code{object$model$curve} is
#'   \code{"exponential"}, \code{"logistic"}, or \code{"richards"},
#'   then \code{top} may also contain \code{"tdoubling"} and \code{"R0"}.
#' @param subset
#'   An expression to be evaluated in the combined model frame
#'   (see \code{\link{egf_combine_frames}}). It must evaluate
#'   to a valid index vector for the rows of the data frame
#'   (see \code{\link{[.data.frame}}), and thus fitting windows.
#'   Confidence intervals are computed only for indexed windows.
#'   The default (\code{\link{NULL}}) is to consider all windows.
#' @param append
#'   An expression indicating variables in the combined model frame
#'   (see \code{\link{egf_combine_frames}}) to be included with the
#'   result. The default (\code{\link{NULL}}) is to append nothing.
#' @param link
#'   A \link{logical} flag. If \code{FALSE}, then confidence intervals
#'   on inverse link-transformed fitted values are computed.
#' @param method
#'   A \link{character} string indicating how confidence intervals
#'   should be computed.
#' @param grid_len (For \code{method = "profile"}.)
#'   A positive integer. Step sizes chosen adaptively by
#'   \code{\link[TMB]{tmbprofile}} will generate approximately
#'   this many points on each side of a profile's minimum point.
#' @param trace (For \code{method != "wald"}.)
#'   A \link{logical} flag.
#'   If \code{TRUE}, then basic tracing messages indicating progress
#'   are printed.
#'   Depending on \code{object$control$trace}, these may be mixed with
#'   optimizer output.
#' @param parallel (For \code{method != "wald"}.)
#'   An \code{"\link{egf_parallel}"} object defining options for \R level
#'   parallelization.
#' @param max_width (For \code{method = "uniroot"}.)
#'   A positive number. \code{\link[TMB]{tmbroot}} will search for roots
#'   in the interval from \code{x-max_width} to \code{x+max_width}, where
#'   \code{x} is the fitted value (link scale).
#' @param breaks,probs (For \code{top = "R0"}.)
#'   Arguments to \code{\link{compute_R0}}.
#' @param .subset,.append
#'   Index vectors to be used (if non-\code{\link{NULL}}) in place of
#'   the result of evaluating \code{subset} and \code{append}.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' Three methods are provided for calculating confidence intervals:
#' \describe{
#' \item{\code{wald}}{
#'   See \code{\link{confint.egf_fitted}}.
#' }
#' \item{\code{profile}}{
#'   See \code{\link{confint.egf_profile}}.
#' }
#' \item{\code{uniroot}}{
#'   Similar to \code{"profile"}, except that the two solutions
#'   of \code{deviance(value) = \link{qchisq}(level, df = 1)} are
#'   approximated by root-finding using \code{\link[TMB]{tmbroot}}
#'   (\code{\link{uniroot}} internally).
#' }
#' }
#' For top level parameters following random effects models,
#' \code{"wald"} returns confidence intervals on individual fitted values,
#' while \code{"profile"} and \code{"uniroot"} return confidence intervals
#' on population fitted values, which are the fixed effects components of
#' individual fitted values.
#'
#' \code{"wald"} assumes e.g., asymptotic normality of the maximum likelihood
#' estimator. \code{"profile"} and \code{"uniroot"} avoid these issues but are
#' typically more expensive, requiring estimation of many restricted models.
#' They are parallelized at the C++ level when there is OpenMP support and
#' \code{object$control$omp_num_threads} is set to an integer greater than 1.
#' If there is no OpenMP support, then computation can still be parallelized
#' at the \R level with appropriate setting of \code{parallel}.
#'
#' See topic \code{\link{egf_eval}} for details on nonstandard evaluation
#' of \code{subset} and \code{append}.
#'
#' @return
#' A \link[=data.frame]{data frame} inheriting from \link{class}
#' \code{"egf_confint"}, with variables:
#' \item{top}{
#'   Top level nonlinear model parameter,
#'   from \code{\link{egf_get_names_top}(object, link = TRUE)}.
#' }
#' \item{ts}{
#'   Time series, from \code{\link{levels}(object$frame_windows$ts)}.
#' }
#' \item{window}{
#'   Fitting window, from \code{\link{levels}(object$frame_windows$window)}.
#' }
#' \item{estimate, lower, upper}{
#'   Fitted value and approximate lower and upper confidence limits.
#' }
#' \code{level} and \code{object$frame_windows} are retained
#' as \link{attributes}.
#'
#' @examples
#' example("egf", "epigrowthfit")
#' zz <- confint(object)
#' str(zz)
#'
#' @seealso \code{\link{plot.egf_confint}}
#' @export
#' @importFrom Matrix KhatriRao
#' @importFrom methods as
#' @importFrom stats qchisq fitted profile confint
#' @import parallel
confint.egf <- function(object,
                        parm,
                        level = 0.95,
                        top = egf_get_names_top(object, link = TRUE),
                        subset = NULL,
                        append = NULL,
                        link = TRUE,
                        method = c("wald", "profile", "uniroot"),
                        grid_len = 12,
                        trace = TRUE,
                        parallel = egf_parallel(),
                        max_width = 7,
                        breaks = NULL,
                        probs = NULL,
                        .subset = NULL,
                        .append = NULL,
                        ...) {
  stop_if_not_number_in_interval(level, 0, 1, "()")
  stop_if_not_true_false(link)
  method <- match.arg(method)
  elu <- c("estimate", "lower", "upper")

  names_top <- names_top_bak <- egf_get_names_top(object, link = TRUE)
  spec <- c("tdoubling", "R0")
  if (object$model$curve %in% c("exponential", "logistic", "richards")) {
    names_top <- c(names_top, spec)
  }

  top <- top_bak <- unique(match.arg(top, names_top, several.ok = TRUE))
  if ("R0" %in% top) {
    stopifnot(
      !is.null(breaks),
      !is.null(probs)
    )
  }
  top[top %in% spec] <- "log(r)"
  top <- unique(top)

  combined <- egf_combine_frames(object)
  subset <- if (is.null(.subset)) substitute(subset) else .subset
  subset <- egf_eval_subset(subset, combined, parent.frame())
  append <- if (is.null(.append)) substitute(append) else .append
  append <- egf_eval_append(append, combined, baseenv())

  if (method == "wald") {
    fo <- fitted(object, top = top, se = TRUE, .subset = subset, .append = append)
    res <- confint(fo, level = level, link = link)

  } else if (method == "profile") {
    po <- profile(object, top = top, .subset = subset, .append = append,
      level_max = level + min(0.01, 0.1 * (1 - level)),
      grid_len = grid_len,
      trace = trace,
      parallel = parallel
    )
    res <- confint(po, level = level, link = link)
    res$linear_combination <- NULL

  } else { # "uniroot"
    stop_if_not_number(max_width, "positive")
    stop_if_not_true_false(trace)
    stopifnot(inherits(parallel, "egf_parallel"))
    n <- length(object$nonrandom)

    p <- length(top)
    N <- length(subset)
    m <- p * N

    f <- factor(object$info$X$top, levels = top)
    J <- as(f, "sparseMatrix")
    X <- object$tmb_out$env$data$X[subset, , drop = FALSE]
    A <- KhatriRao(J, X)
    if (egf_has_random(object)) {
      index <- grepl("^beta\\[", names(object$best)[object$nonrandom])
      A@p <- c(0L, replace(rep_len(NA_integer_, n), index, A@p[-1L]))
      A@p <- locf(A@p)
      A@Dim <- c(m, n)
    }

    a <- lapply(seq_len(m), function(i) A[i, ])
    target <- 0.5 * qchisq(level, df = 1) # y := diff(nll) = 0.5 * deviance
    sd.range <- max_width
    nomp <- object$control$omp_num_threads

    do_uniroot <- function(i, a) {
      if (trace) {
        cat(sprintf("Computing confidence interval %d of %d...\n", i, m))
      }
      TMB::tmbroot(obj, lincomb = a, target = target, sd.range = sd.range, trace = FALSE)
    }

    if (parallel$method == "snow") {
      environment(do_uniroot) <- .GlobalEnv

      ## Reconstruct list of arguments to 'MakeADFun' from object internals
      ## for retaping
      args <- egf_tmb_remake_args(object)

      ## Retrieve path to shared object for loading
      dll <- system.file("libs", TMB::dynlib("epigrowthfit"), package = "epigrowthfit", mustWork = TRUE)

      if (is.null(parallel$cl)) {
        cl <- do.call(makePSOCKcluster, parallel$args)
        on.exit(stopCluster(cl), add = TRUE)
      } else {
        cl <- parallel$cl
      }
      clusterExport(cl,
        varlist = c("nomp", "args", "trace", "m", "target", "sd.range"),
        envir = environment()
      )
      clusterEvalQ(cl, {
        dyn.load(dll)
        if (TMB::openmp(n = NULL) > 0L) {
          TMB::openmp(n = nomp)
        }
        obj <- do.call(TMB::MakeADFun, args)
      })
      res <- clusterMap(cl, do_uniroot, i = seq_len(m), a = a)
    } else {
      obj <- object$tmb_out
      if ((onomp <- TMB::openmp(n = NULL)) > 0L) {
        TMB::openmp(n = nomp)
        on.exit(TMB::openmp(n = onomp), add = TRUE)
      }
      if (given_outfile <- nzchar(parallel$outfile)) {
        outfile <- file(parallel$outfile, open = "wt")
        sink(outfile, type = "output")
        sink(outfile, type = "message")
      }
      res <- tryCatch(
        switch(parallel$method,
          multicore = do.call(mcMap, c(list(f = do_uniroot, i = seq_len(m), a = a), parallel$args)),
          serial = Map(do_uniroot, i = seq_len(m), a = a)
        ),
        error = function(cond) {
          if (given_outfile) {
            sink(type = "message")
            sink(type = "output")
          }
          stop(cond)
        }
      )
    }

    res <- data.frame(
      top = rep(factor(top, levels = names_top), each = N),
      ts = object$frame_windows$ts[subset],
      window = object$frame_windows$window[subset],
      estimate = as.numeric(A %*% object$best[object$nonrandom]),
      `colnames<-`(do.call(rbind, res), c("lower", "upper")),
      combined[subset, append, drop = FALSE],
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    if (!link) {
      res[elu] <- in_place_ragged_apply(res[elu], res$top,
        f = lapply(egf_link_extract(levels(res$top)), egf_link_match, inverse = TRUE)
      )
      levels(res$top) <- egf_link_remove(levels(res$top))
    }
  }

  if (any(top_bak %in% spec)) {
    s <- if (link) "log(r)" else "r"
    res_r <- res[res$top == s, , drop = FALSE]
    if (link) {
      res_r[elu] <- exp(res_r[elu])
    }

    if ("tdoubling" %in% top_bak) {
      eul <- elu[c(1L, 3L, 2L)]
      res_tdoubling <- res_r
      res_tdoubling[elu] <- log(2) / res_r[eul]
      res_tdoubling$top <- factor("tdoubling")
      res <- rbind(res, res_tdoubling)
    }

    if ("R0" %in% top_bak) {
      res_R0 <- res_r
      res_R0[elu] <- lapply(res_r[elu], compute_R0, breaks = breaks, probs = probs)
      res_R0$top <- factor("R0")
      res <- rbind(res, res_R0)
    }

    if (!"log(r)" %in% top_bak) {
      res <- res[res$top != s, , drop = FALSE]
    }
  }

  row.names(res) <- NULL
  attr(res, "level") <- level
  attr(res, "frame_windows") <- object$frame_windows # for 'plot.egf_confint'
  class(res) <- c("egf_confint", "data.frame")
  res
}

#' Plot confidence intervals
#'
#' A method for graphically comparing confidence intervals
#' on fitted values of top level nonlinear model parameters.
#'
#' @param x
#'   An \code{"\link[=confint.egf]{egf_confint}"} object.
#' @param type
#'   A \link{character} string determining how confidence intervals
#'   are displayed (see Details).
#' @param time_as
#'   A \link{character} string indicating how time is displayed
#'   on the bottom axis. The options are: as is (\code{"numeric"})
#'   and with a calendar (\code{"Date"}). In the latter case,
#'   numeric times are interpreted as numbers of days since
#'   \code{1970-01-01 00:00:00}. [\code{type = "boxes"} only.]
#' @param per_plot
#'   A positive integer. One plot will display at most this many
#'   confidence intervals or time series (depending on \code{type}).
#' @param subset
#'   An expression to be evaluated in \code{x}.
#'   It must evaluate to a valid index vector for the rows of \code{x}
#'   (see \code{\link{[.data.frame}}).
#'   Only indexed confidence intervals are plotted.
#'   The default (`NULL`) is to plot all confidence intervals.
#' @param order
#'   An expression to be evaluated in \code{x},
#'   typically a call to \code{\link{order}},
#'   determining the order in which confidence intervals
#'   or time series (depending on \code{type}) are plotted.
#'   It must evaluate to a permutation of \code{\link{seq_len}(\link{nrow}(x))}.
#'   The default (\code{\link{NULL}}) is equivalent to \code{seq_len(nrow(x))}.
#' @param label
#'   An expression to be evaluated in \code{x}, typically a \link{factor},
#'   \link{interaction} of factors, or \link{call} to \code{\link{sprintf}}.
#'   It is used to create appropriate \eqn{y}-axis labels for confidence
#'   intervals when \code{type = "bars"} and appropriate panel labels for
#'   time series when \code{type = "boxes"}. The default (\code{\link{NULL}})
#'   is to take labels from \code{x$window} and \code{x$ts}, respectively.
#' @param main
#'   An expression or character string indicating a plot title,
#'   to be recycled for all plots.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' \code{type = "bars"} creates a one-dimensional plot with
#' confidence intervals drawn as stacked horizontal line segments.
#'
#' \code{type = "boxes"} creates a two-dimensional plot for each
#' time series (level of \code{x$ts}), each containing shaded boxes.
#' Projection of boxes onto the horizontal and vertical axes yields
#' fitting windows and corresponding confidence intervals, respectively.
#'
#' If an endpoint of a confidence interval is \code{\link{NaN}},
#' then dashed lines are drawn from the point estimate to
#' the boundary of the plotting region to indicate missingness.
#'
#' See topic \code{\link{egf_eval}} for details on nonstandard evaluation
#' of \code{subset}, \code{order}, and \code{label}.
#'
#' @return
#' \code{\link{NULL}} (invisibly).
#'
#' @examples
#' example("confint.egf", "epigrowthfit")
#'
#' op <- par(mar = c(4.5, 4, 2, 1), oma = c(0, 0, 0, 0))
#' plot(zz, type = "bars")
#' par(op)
#'
#' op <- par(mar = c(0.2, 0, 0.2, 0), oma = c(4.5, 6, 2, 1), las = 1)
#' plot(zz, type = "boxes")
#' par(op)
#'
#' @export
plot.egf_confint <- function(x,
                             type = c("bars", "boxes"),
                             time_as = c("Date", "numeric"),
                             per_plot = switch(type, bars = 12L, boxes = 3L),
                             subset = NULL,
                             order = NULL,
                             main = NULL,
                             label = NULL,
                             ...) {
  type <- match.arg(type)
  time_as <- match.arg(time_as)
  stop_if_not_integer(per_plot, "positive")

  subset <- egf_eval_subset(substitute(subset), x, parent.frame())
  if (length(subset) == 0L) {
    stop("'subset' indexes zero confidence intervals, so there is nothing to plot.")
  }
  label <- egf_eval_label(substitute(label), x, parent.frame())
  order <- egf_eval_order(substitute(order), x, parent.frame())
  subset <- order[order %in% subset]

  a <- attributes(x)
  if (is.null(main)) {
    s <- switch(type, bars = "fitting window", boxes = "time series")
    main <- sprintf("%.3g%% confidence intervals by %s", 100 * a$level, s)
  }
  if (is.null(label)) {
    s <- switch(type, bars = "window", boxes = "ts")
    label <- as.character(x[[s]])
  }

  nx <- c("top", "ts", "window", "estimate", "lower", "upper")
  x <- x[nx]
  x$label <- label
  x <- x[subset, , drop = FALSE]
  x$top <- factor(x$top)

  if (type == "bars") {
    plot.egf_confint.bars(x, per_plot = per_plot, main = main)
  } else {
    i <- match(x$window, a$frame_windows$window, 0L)
    frame_windows <- a$frame_windows[i, c("start", "end"), drop = FALSE]
    x <- data.frame(x, frame_windows)
    plot.egf_confint.boxes(x, time_as = time_as, per_plot = per_plot, main = main)
  }
}

#' @import graphics
plot.egf_confint.bars <- function(x, per_plot, main) {
  gp <- par()
  sfcex <- get_space_filling_cex(x$label, target = 0.95 * max(0, gp$mar[2L] - 0.25), units = "lines")

  for (top in levels(x$top)) { # loop over parameters
    data <- x[x$top == top, , drop = FALSE]
    argna <- lapply(data[c("lower", "upper")], is.na)
    n <- nrow(data)
    xlab <- top
    xlim <- range(data[c("estimate", "lower", "upper")], na.rm = TRUE)

    i <- 0L
    while (i < n) { # loop over plots
      k <- i + seq_len(min(per_plot, n - i))
      plot.new()
      plot.window(xlim = xlim, ylim = c(per_plot + 1, 0), xaxs = "r", yaxs = "i")
      gp$usr <- par("usr")
      abline(v = axTicks(side = 1), lty = 3, col = "grey75")
      segments(
        x0  = replace(data$lower[k], argna$lower[k], gp$usr[1L]),
        x1  = data$estimate[k],
        y0  = seq_along(k),
        y1  = seq_along(k),
        lty = 1 + as.numeric(argna$lower[k]),
        lwd = 2 - as.numeric(argna$lower[k])
      )
      segments(
        x0  = data$estimate[k],
        x1  = replace(data$upper[k], argna$upper[k], gp$usr[2L]),
        y0  = seq_along(k),
        y1  = seq_along(k),
        lty = 1 + as.numeric(argna$upper[k]),
        lwd = 2 - as.numeric(argna$upper[k])
      )
      points(
        x   = data$estimate[k],
        y   = seq_along(k),
        pch = 21,
        bg  = "grey80"
      )
      box()
      axis(side = 1)
      mtext(
        text = data$label[k],
        side = 2,
        line = 0.25,
        at = seq_along(k),
        las = 1,
        adj = 1,
        padj = 0.5,
        cex = gp$cex * sfcex
      )
      title(xlab = xlab)
      title(main, adj = 0)
      i <- i + per_plot
    } # loop over plots
  } # loop over parameters

  invisible(NULL)
}

#' @import graphics
plot.egf_confint.boxes <- function(x, time_as, per_plot, main) {
  op <- par(mfrow = c(per_plot, 1))
  on.exit(par(op))
  gp <- par()
  xlim <- range(x[c("start", "end")])

  for (top in levels(x$top)) { # loop over parameters
    data <- x[x$top == top, , drop = FALSE]
    ylim <- range(data[c("estimate", "lower", "upper")], na.rm = TRUE)
    ylab <- top

    data <- split(data, factor(data$ts, levels = as.character(unique(data$ts))))
    n <- length(data)

    i <- 0L
    while (i < n) { # loop over pages
      for (k in i + seq_len(min(per_plot, n - i))) { # loop over plots
        plot.new()
        plot.window(xlim = xlim, ylim = ylim)

        v <- Daxis(side = 1, major = NULL, minor = NULL)$minor
        abline(v = v, lty = 3, col = "grey75")

        gp$usr <- par("usr")
        argna <- lapply(data[[k]][c("lower", "upper")], is.na)
        data[[k]]$lower[argna$lower] <- gp$usr[3L]
        data[[k]]$upper[argna$upper] <- gp$usr[4L]

        for (l in seq_len(nrow(data[[k]]))) { # loop over confidence intervals
          rect(
            xleft   = data[[k]]$start[l],
            xright  = data[[k]]$end[l],
            ybottom = data[[k]]$lower[l],
            ytop    = data[[k]]$estimate[l],
            col     = if (argna$lower[l]) NA else "grey80",
            border  = "grey50",
            lty     = if (argna$lower[l]) 2 else 1
          )
          rect(
            xleft   = data[[k]]$start[l],
            xright  = data[[k]]$end[l],
            ybottom = data[[k]]$estimate[l],
            ytop    = data[[k]]$upper[l],
            col     = if (argna$upper[l]) NA else "grey80",
            border  = "grey50",
            lty     = if (argna$upper[l]) 2 else 1
          )
          segments(
            x0  = data[[k]]$start[l],
            x1  = data[[k]]$end[l],
            y0  = data[[k]]$estimate[l],
            y1  = data[[k]]$estimate[l],
            col = "grey50",
            lwd = 2
          )
        } # loop over confidence intervals

        p <- diff(grconvertY(c(0, 0.08), "npc", "inches"))
        px <- diff(grconvertX(c(0, p), "inches", "user"))
        py <- diff(grconvertY(c(0, p), "inches", "user"))

        ocex <- par(cex = 1)
        text(
          x = gp$usr[1L] + px,
          y = gp$usr[4L] - py,
          labels = data[[k]]$label[1L],
          adj = c(0, 1),
          cex = gp$cex.lab
        )
        par(ocex)

        box()
        axis(side = 2)
      } # loop over plots

      ocex <- par(cex = 1)
      if (time_as == "numeric") {
        axis(side = 1)
      } else {
        Daxis(
          side = 1,
          major = list(cex.axis = 1.2 * gp$cex.axis, mgp = gp$mgp + c(0, 1.5, 0), tick = FALSE)
        )
      }
      par(ocex)

      ## Hack to avoid dealing with normalized figure coordinates
      par(new = TRUE, mfrow = c(1, 1), mai = gp$mai + gp$omi, omi = c(0, 0, 0, 0))
      plot.window(xlim = xlim, ylim = c(0, 1))
      title(main = main, adj = 0)
      title(ylab = ylab, adj = 0.5)
      par(new = FALSE, mfrow = gp$mfrow, mai = gp$mai, omi = gp$omi)

      i <- i + per_plot
    } # loop over pages
  } # loop over parameters

  invisible(NULL)
}
