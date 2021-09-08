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
#' @param parallel (For \code{method = "profile"}.)
#'   An \code{"\link{egf_parallel}"} object defining parallelization
#'   options.
#' @param trace (For \code{method != "wald"}.)
#'   A \link{logical} flag.
#'   If \code{TRUE}, then basic tracing messages indicating progress
#'   are printed.
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
#' whereas \code{"profile"} and \code{"uniroot"} return confidence intervals
#' population fitted values, which are the fixed effects components
#' of individual fitted values.
#'
#' \code{"wald"} requires minimal computation time but assumes,
#' e.g., asymptotic normality of the maximum likelihood estimator.
#' A further limitation of \code{"wald"} is functional non-invariance.
#' \code{"profile"} and \code{"uniroot"} avoid these issues but are
#' much slower, requiring estimation of restricted models. Of the two,
#' \code{"profile"} is more robust.
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
#' @seealso \code{plot.egf_confint}
#' @export
#' @importFrom Matrix KhatriRao
#' @importFrom methods as
#' @importFrom stats qchisq fitted profile confint
#' @importFrom TMB tmbroot openmp
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
                        parallel = egf_parallel(),
                        trace = TRUE,
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
    ft <- fitted(object, top = top, se = TRUE, .subset = subset, .append = append)
    res <- confint(ft, level = level, link = link)

  } else if (method == "profile") {
    pf <- profile(object, top = top, .subset = subset, .append = append,
      max_level = level + min(0.01, 0.1 * (1 - level)),
      grid_len = grid_len,
      parallel = parallel,
      trace = trace
    )
    res <- confint(pf, level = level, link = link)
    res$linear_combination <- NULL

  } else { "uniroot"
    stop_if_not_number(max_width, "positive")

    p <- length(top)
    N <- sum(subset)
    m <- p * N
    n <- length(object$nonrandom)

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
    A_rows <- lapply(seq_len(m), function(i) A[i, ])

    tmbroot_args <- list(
      obj = object$tmb_out,
      target = qchisq(level, df = 1) / 2, # y := diff(nll) = deviance / 2
      sd.range = max_width,
      trace = FALSE
    )
    omp_num_threads <- object$control$omp_num_threads

    do_uniroot <- function(r, i) {
      if (trace) {
        cat(sprintf("Computing confidence interval %d of %d...\n", i, m))
      }
      on <- openmp(n = NULL)
      if (on > 0L) {
        openmp(n = omp_num_threads)
        on.exit(openmp(n = on))
      }
      tmbroot_args$lincomb <- r
      do.call(tmbroot, tmbroot_args)
    }

    if (parallel == "snow") {
      environment(do_uniroot) <- .GlobalEnv # see comment in R/boot.R
      if (is.null(parallel$cl)) {
        cl <- do.call(makePSOCKcluster, parallel$options)
        on.exit(stopCluster(cl))
      } else {
        cl <- parallel$cl
      }
      clusterEvalQ(cl, library("TMB"))
      clusterExport(cl,
        varlist = c("tmbroot_args", "trace", "m", "omp_num_threads"),
        envir = environment()
      )
      vl <- clusterMap(cl, do_uniroot, r = A_rows, i = seq_len(m))

    } else {
      if (nzchar(parallel$outfile)) {
        outfile <- file(parallel$outfile, open = "wt")
        sink(outfile, type = "output")
        sink(outfile, type = "message")
      }
      vl <- switch(parallel$method,
        multicore = do.call(mcMap, c(list(f = do_uniroot, r = A_rows, i = seq_len(m)), parallel$options)),
        serial = Map(do_uniroot, r = A_rows, i = seq_len(m))
      )
      if (nzchar(parallel$outfile)) {
        sink(type = "output")
        sink(type = "message")
      }
    }

    res <- data.frame(
      top = rep(factor(top, levels = names_top), each = N),
      ts = object$frame_windows$ts[subset],
      window = object$frame_windows$window[subset],
      estimate = as.numeric(A %*% object$best[object$nonrandom]),
      do.call(rbind, vl), # "lower" "upper"
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

# #' Plot confidence intervals
# #'
# #' A method for graphically comparing confidence intervals
# #' on fitted values of nonlinear model parameters across
# #' fitting windows.
# #'
# #' @param x
# #'   An `"egf_confint"` object returned by [confint.egf()].
# #' @param type
# #'   A character string determining how confidence intervals
# #'   are displayed (see Details).
# #' @param subset
# #'   An expression to be evaluated in `x`. Must evaluate to
# #'   a logical vector indexing rows of `x`. Only indexed
# #'   confidence intervals are plotted. The default (`NULL`)
# #'   is to plot all confidence intervals.
# #' @param order
# #'   An expression to be evaluated in `x`, typically a call
# #'   to [order()], determining the order in which confidence
# #'   intervals or time series (depending on `type`) are plotted.
# #'   Must evaluate to a permutation of `seq_len(nrow(x))`.
# #'   The default (`NULL`) is equivalent to `seq_len(nrow(x))`.
# #' @param per_plot
# #'   A positive integer. One plot will display at most this many
# #'   confidence intervals or time series (depending on `type`).
# #' @param main
# #'   An expression or character string indicating a plot title,
# #'   to be recycled for all plots.
# #' @param label
# #'   An expression to be evaluated in `x`, typically a factor,
# #'   interaction of factor, or call to [sprintf()]. Used to
# #'   create appropriate _y_-axis labels for confidence intervals
# #'   when `type = "bars"` and approprate panel labels for
# #'   time series when `type = "boxes"`. The default (`NULL`)
# #'   is to take labels from `x$window` and `x$ts`, respectively.
# #' @param ...
# #'   Unused optional arguments.
# #'
# #' @details
# #' `type = "bars"` creates a one-dimensional plot with
# #' confidence intervals drawn as stacked horizontal line segments.
# #'
# #' `type = "boxes"` creates a two-dimensional plot for each
# #' time series (level of `x$ts`), each containing shaded boxes.
# #' Projection of boxes onto the horizontal and vertical axes
# #' yields fitting windows and corresponding confidence intervals,
# #' respectively.
# #'
# #' If an endpoint of a confidence interval is `NA`, then dashed
# #' lines are drawn from the point estimate to the boundary of the
# #' plotting region to indicate missingness.
# #'
# #' @return
# #' `NULL` (invisibly).
# #'
# #' @export
# plot.egf_confint <- function(x,
#                              type = c("bars", "boxes"),
#                              subset = NULL,
#                              order = NULL,
#                              per_plot = switch(type, bars = 12L, 4L),
#                              main = NULL,
#                              label = NULL,
#                              ...) {
#   type <- match.arg(type)
#   stop_if_not_integer(per_plot, kind = "positive")
#
#   subset <- subset_to_index(substitute(subset), x, parent.frame())
#   order <- order_to_index(substitute(order), x, parent.frame())
#
#   a <- attributes(x)
#   x <- x[order[order %in% subset], , drop = FALSE]
#
#   if (is.null(main)) {
#     s <- switch(type, bars = "fitting window", "time series")
#     main <- sprintf("%.3g%% CI by %s", 100 * a$level, s)
#   }
#
#   label <- label_to_character(substitute(label), x, parent.frame())
#   if (is.null(label)) {
#     s <- switch(type, bars = "window", "ts")
#     label <- as.character(x[[s]])
#   }
#
#   nx <- c("par", "ts", "window", "estimate", "lower", "upper")
#   x <- data.frame(x[nx], label, stringsAsFactors = FALSE)
#
#   if (type == "bars") {
#     do_bars_plot(x, per_plot = per_plot, main = main)
#   } else {
#     shift <- min(a$endpoints$start)
#     origin <- attr(a$endpoints, "origin") + shift
#     iep <- match(x$window, a$endpoints$window, 0L)
#     endpoints <- a$endpoints[iep, c("start", "end"), drop = FALSE] - shift
#     x <- data.frame(x, endpoints)
#
#     do_boxes_plot(x, origin = origin, per_plot = per_plot, main = main)
#   }
#   invisible(NULL)
# }
#
# #' @keywords internal
# #' @import graphics
# do_bars_plot <- function(x, per_plot, main) {
#   mar <- c(3.5, 5, 1.5, 1)
#   csi <- par("csi")
#   op <- par(mar = mar)
#   on.exit(par(op))
#
#   x_split <- split(x, factor(x$par, levels = unique(x$par)))
#   nxs <- names(x_split)
#
#   yax_cex <- get_yax_cex(x$label, mex = 0.92 * mar[2L], cex = 0.8, font = 1, csi = csi)
#   yax_cex <- min(0.8, yax_cex)
#
#   for (i in seq_along(x_split)) { # loop over nonlinear model parameters
#     xi <- x_split[[i]]
#     xlab <- nxs[i]
#     xlim <- range(xi[c("estimate", "lower", "upper")], na.rm = TRUE)
#     lna <- is.na(xi$lower)
#     una <- is.na(xi$upper)
#
#     K <- 0L
#     while (K < nrow(xi)) { # loop over plots
#       k <- K + seq_len(min(per_plot, nrow(xi) - K))
#       plot.new()
#       plot.window(xlim = xlim, ylim = c(per_plot + 1, 0), xaxs = "r", yaxs = "i")
#       usr <- par("usr")
#       abline(v = axTicks(side = 1), lty = 3, col = "grey75")
#       segments(
#         x0  = replace(xi$lower[k], lna[k], usr[1L]),
#         x1  = xi$estimate[k],
#         y0  = seq_along(k),
#         y1  = seq_along(k),
#         lty = c(1, 2)[1L + lna[k]],
#         lwd = c(2, 1)[1L + lna[k]]
#       )
#       segments(
#         x0  = xi$estimate[k],
#         x1  = replace(xi$upper[k], una[k], usr[2L]),
#         y0  = seq_along(k),
#         y1  = seq_along(k),
#         lty = c(1, 2)[1L + una[k]],
#         lwd = c(2, 1)[1L + una[k]]
#       )
#       points(
#         x   = xi$estimate[k],
#         y   = seq_along(k),
#         pch = 21,
#         bg  = "grey75"
#       )
#       box()
#       axis(
#         side = 1,
#         tcl = -0.4,
#         mgp = c(3, 0.3, 0),
#         cex.axis = yax_cex
#       )
#       axis(
#         side = 2,
#         at = seq_along(k),
#         labels = xi$label[k],
#         tick = FALSE,
#         las = 1,
#         mgp = c(3, 0.2, 0),
#         cex.axis = 0.8
#       )
#       title(xlab = xlab, line = 2, cex.lab = 0.9)
#       title(main, line = 0.25, adj = 0, cex.main = 0.9)
#       K <- K + per_plot
#     } # loop over plots
#   } # loop over nonlinear model parameters
#
#   invisible(NULL)
# }
#
# #' @keywords internal
# #' @import graphics
# do_boxes_plot <- function(x, origin, per_plot, main) {
#   op <- par(
#     mfrow = c(per_plot, 1),
#     mar = c(0, 3, 0.25, 0.5),
#     oma = c(2.5, 1.5, 1.5, 0)
#   )
#   on.exit(par(op))
#
#   x_split <- split(x, factor(x$par, levels = unique(x$par)))
#   nxs <- names(x_split)
#
#   xlim <- c(0, max(x$end))
#   xax_at <- daxis(origin = origin + 1, plot = FALSE)
#
#   for (i in seq_along(x_split)) { # loop over nonlinear model parameters
#     xi <- x_split[[i]]
#     xi_split <- split(xi, factor(xi$ts, levels = unique(xi$ts)))
#     ylab <- nxs[i]
#     ylim <- range(xi[c("estimate", "lower", "upper")], na.rm = TRUE)
#
#     K <- 0L
#     while (K < length(xi_split)) { # loop over plots
#       for (k in K + seq_len(min(per_plot, length(xi_split) - K))) { # loop over panels
#         xik <- xi_split[[k]]
#
#         plot.new()
#         plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "r")
#         gp <- par(c("usr", "cex", "mai", "omi", "pin", "din"))
#         abline(v = xax_at, lty = 3, col = "grey75")
#
#         xik$lower[lna <- is.na(xik$lower)] <- gp$usr[3L]
#         xik$upper[una <- is.na(xik$upper)] <- gp$usr[4L]
#
#         for (l in seq_len(nrow(xik))) { # loop over boxes
#           rect(
#             xleft   = xik$start[l],
#             xright  = xik$end[l],
#             ybottom = xik$lower[l],
#             ytop    = xik$estimate[l],
#             col     = if (lna[l]) NA else "grey75",
#             border  = "grey50",
#             lty     = if (lna[l]) 2 else 1
#           )
#           rect(
#             xleft   = xik$start[l],
#             xright  = xik$end[l],
#             ybottom = xik$estimate[l],
#             ytop    = xik$upper[l],
#             col     = if (una[l]) NA else "grey75",
#             border  = "grey50",
#             lty     = if (una[l]) 2 else 1
#           )
#           segments(
#             x0  = xik$start[l],
#             x1  = xik$end[l],
#             y0  = xik$estimate[l],
#             y1  = xik$estimate[l],
#             col = "grey50",
#             lwd = 2
#           )
#         } # loop over boxes
#
#         box()
#         axis(
#           side = 2,
#           las = 1,
#           mgp = c(3, 0.7, 0),
#           cex.axis = 0.8
#         )
#         text(
#           x = gp$usr[1L] + 0.075 * (gp$usr[2L] - gp$usr[1L]) * (gp$pin[2L] / gp$pin[1L]),
#           y = gp$usr[4L] - 0.075 * (gp$usr[4L] - gp$usr[3L]),
#           labels = xik$label[1L],
#           adj = c(0, 1),
#           cex = 0.8
#         )
#       } # loop over panels
#
#       mtext(main,
#         side = 3,
#         line = 0,
#         outer = TRUE,
#         at = gp$mai[2L] / (gp$din[1L] - sum(gp$omi[c(2L, 4L)])),
#         adj = 0,
#         cex = 0.9,
#         font = 2
#       )
#       mtext(ylab,
#         side = 2,
#         line = 0,
#         outer = TRUE,
#         at = 0.5,
#         adj = 0.5,
#         cex = 0.9
#       )
#       daxis(
#         origin = origin + 1,
#         minor = list(
#           mgp = c(3, 0.25, 0), xpd = TRUE,
#           lwd = 0, lwd.ticks = 1, tcl = -0.2,
#           gap.axis = 0, cex.axis = 0.8 / gp$cex
#         ),
#         major = list(
#           mgp = c(3, 1.25, 0), xpd = TRUE,
#           lwd = 0, lwd.ticks = 0, tcl = 0,
#           gap.axis = 0, cex.axis = 0.9 / gp$cex
#         )
#       )
#       K <- K + per_plot
#     } # loop over plots
#   } # loop over nonlinear model parameters
#
#   invisible(NULL)
# }
