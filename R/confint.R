#' Confidence intervals on fitted values
#'
#' Computes confidence intervals on fitted values of response variables
#' (see [fitted.egf()]) and, where appropriate, the basic reproduction
#' number and initial doubling time in days.
#'
#' @param object
#'   An `"egf"` object returned by [egf()].
#' @param parm
#'   A character vector giving a subset of the response variables
#'   listed in `get_par_names(fitted, link = TRUE)`.
#'   If `object$curve` is an element of
#'   `c("exponential", "logistic", "richards")`,
#'   then `parm` may also contain `"R0"` and `"tdoubling"`.
#' @param level
#'   A number in the interval (0,1). The desired confidence level.
#' @param link
#'   A logical scalar. If `FALSE`, then confidence intervals on
#'   inverse-link transformed fitted values are returned.
#' @param method
#'   A character string indicating how confidence intervals should
#'   be calculated (see Details).
#' @param max_width
#'   A positive numeric vector, to be recycled to length `length(parm)`.
#'   [TMB::tmbroot()] will search for roots in the interval from
#'   `x-max_width` to `x+max_width`, where `x` is the fitted value
#'   (link scale) of the response variable of interest.
#' @param breaks,probs
#'   Arguments to `compute_R0()`, necessary if `"R0 %in% parm"`.
#' @param ...
#'   Unused optional arguments.
#' @inheritParams profile.egf
#'
#' @details
#' Three methods are provided for calculating confidence intervals:
#' \describe{
#' \item{`wald`}{
#'   Confidence limits on a fitted value (link scale) are computed as
#'   `estimate + c(-1, 1) * sqrt(q) * se`,
#'   where `q = qchisq(level, df = 1)` and `estimate` and `se` are the
#'   fitted value and its approximate (delta method) standard error.
#' }
#' \item{`profile`}{
#'   A discrete approximation of the negative log likelihood profile
#'   is computed using [TMB::tmbprofile()] then linearly interpolated
#'   using [TMB::confint.tmbprofile()] to approximate the two solutions
#'   of `deviance(x) = qchisq(level, df = 1)`. These are the lower
#'   and upper confidence limits on the fitted value (link scale).
#'   See Wilks' theorem.
#' }
#' \item{`uniroot`}{
#'   Similar to `"profile"`, except that the two solutions of
#'   `deviance(x) = qchisq(level, df = 1)` are approximated
#'   by root-finding using [TMB::tmbroot()]. [stats::uniroot()]
#'   is called internally.
#' }
#' }
#' `"wald"` is the only method available for response variables
#' following a mixed (rather than fixed) effects model.
#'
#' `"wald"` requires minimal computation time but assumes, e.g.,
#' asymptotic normality of the maximum likelihood estimator.
#' A further limitation of `"wald"` is functional non-invariance.
#' `"profile"` and `"uniroot"` avoid these issues but are much
#' slower, requiring estimation of restricted models. Of the two,
#' `"uniroot"` is typically less robust.
#'
#' @return
#' A data frame inheriting from class `"egf_confint"`,
#' with one row per element of `parm` _per fitting window_.
#'
#' The first variable, `par`, is a factor whose levels are
#' `parm` (with prefixes stripped if `link = FALSE`). The next
#' `length(object$frame)-2` variables specify levels for each
#' factor in `object$frame[-(1:2)]` (and thus fitting windows).
#' The remaining variables, `estimate`, `lower`, and `upper`,
#' list fitted values and the corresponding lower and upper
#' confidence limits.
#'
#' @seealso [plot.egf_confint()]
#' @export
#' @importFrom stats qchisq confint profile
#' @import parallel
confint.egf <- function(object, parm = get_par_names(object),
                        level = 0.95, link = TRUE,
                        method = c("wald", "profile", "uniroot"),
                        grid_len = 12, max_width = 7,
                        trace = TRUE,
                        parallel = c("serial", "multicore", "snow"),
                        cores = getOption("egf.cores", 2L),
                        outfile = NULL,
                        cl = NULL,
                        breaks = NULL, probs = NULL, ...) {
  pn <- pn0 <- get_par_names(object, link = TRUE)
  if (object$curve %in% c("exponential", "logistic", "richards")) {
    pn0 <- c("R0", "tdoubling", pn0)
  }
  stop_if_not(
    is.character(parm),
    length(parm) > 0L,
    parm %in% c(pn0, remove_link_string(pn0)),
    m = paste0(
      "`parm` must be a subset of:\n",
      paste(sprintf("\"%s\"", pn0), collapse = ", ")
    )
  )
  if ("R0" %in% parm && (is.null(breaks) || is.null(probs))) {
    stop("`parm = \"R0\"` requires non-NULL `breaks`, `probs`.\n",
         "See `help(\"compute_R0\")`.")
  }
  stop_if_not_in_0_1(level)
  stop_if_not_tf(link)
  method <- match.arg(method)

  parm <- parm0 <- add_link_string(unique(remove_link_string(parm)))
  parm[parm %in% c("R0", "tdoubling")] <- "log_r"
  parm <- unique(parm)
  fr <- object$frame[!duplicated(object$index), -(1:2), drop = FALSE]

  if (method == "wald") {
    estimate <- object$report$Y_short_as_vector$estimate
    se <- object$report$Y_short_as_vector$se
    dim(estimate) <- dim(se) <- c(nrow(fr), length(pn))

    q <- qchisq(level, df = 1)
    get_elu <- function(j) {
      elu <- estimate[, j] + outer(se[, j], c(0, -1, 1) * sqrt(q))
      colnames(elu) <- c("estimate", "lower", "upper")
      elu
    }
    ml <- lapply(match(parm, pn), get_elu)

    if (!link) {
      delink <- function(x, s) get_inverse_link(s)(x)
      ml <- Map(delink, x = ml, s = extract_link_string(parm))
      parm <- remove_link_string(parm)
    }
    out <- cbind(
      par = rep(factor(parm, levels = parm), each = nrow(fr)),
      do.call(rbind, rep.int(list(fr), length(parm))),
      do.call(rbind, ml)
    )
  } else {
    parallel <- match.arg(parallel)
    check_parallel(parallel, cores, outfile, cl)
    rid <- object$tmb_args$data$rid
    stop_if_not(
      colSums(rid[, parm, drop = FALSE]) == 0L,
      m = paste0(
        "Unable to compute likelihood profiles\n",
        "of sums of fixed and random effects:\n",
        paste(colnames(rid)[colSums(rid) > 0L], collapse = ", ")
      )
    )

    if (method == "profile") {
      p <- profile(object,
        parm = parm,
        decontrast = TRUE,
        max_level = level + min(0.01, 0.1 * (1 - level)),
        grid_len = grid_len,
        trace = trace,
        parallel = parallel,
        cores = cores,
        outfile = outfile,
        cl = cl
      )
      out_point_five <- confint(p, level = level)
    } else if (method == "uniroot") {
      stop_if_not(
        is.numeric(max_width),
        length(max_width) > 0L,
        max_width > 0,
        !is.infinite(max_width),
        m = "`max_width` must be a numeric vector with positive elements."
      )

      lin_comb <- make_lin_comb_for_parm(object, parm) # see R/profile.R
      lcl <- lapply(seq_len(nrow(lin_comb)), function(i) lin_comb[i, ])
      wl <- rep(max_width, length.out = length(parm))
      ytol <- qchisq(level, df = 1) / 2 # y := diff(nll) = deviance / 2

      n <- length(lcl)
      iofn <- function(i) sprintf("%*d of %d", nchar(n), i, n)
      f <- function(lc, w, i) {
        if (trace) cat("Computing confidence interval", iofn(i), "...\n")
        TMB::tmbroot(object$tmb_out, target = ytol, lincomb = lc,
                     sd.range = w, trace = FALSE)
      }

      if (parallel == "snow") {
        environment(f) <- environment(iofn) <- .GlobalEnv # see R/boot.R

        if (is.null(cl)) {
          if (is.null(outfile)) {
            outfile <- ""
          }
          cl <- makePSOCKcluster(cores, outfile = outfile)
          on.exit(stopCluster(cl))
        }
        clusterExport(cl,
          varlist = c("ytol", "iofn", "n"),
          envir = environment()
        )
        cil <- clusterMap(cl, f, lc = lcl, w = wl, i = seq_along(lcl))
      } else {
        if (!is.null(outfile)) {
          sink(outfile, type = "output")
          sink(outfile, type = "message")
        }
        cil <- switch(parallel,
          multicore = mcmapply(f, lc = lcl, w = wl, i = seq_along(lcl),
                               SIMPLIFY = FALSE, mc.cores = cores),
          serial = Map(f, lc = lcl, w = wl, i = seq_along(lcl))
        )
        if (!is.null(outfile)) {
          sink(type = "output")
          sink(type = "message")
        }
      }

      estimate <- as.vector(lin_comb %*% object$par[object$nonrandom])
      ci <- do.call(rbind, cil)
      colnames(ci) <- c("lower", "upper")
      out_point_five <- data.frame(
        name = rownames(lin_comb),
        estimate,
        ci,
        stringsAsFactors = FALSE
      )
    }

    ## To create output matching `method = "wald"`, we need to map
    ## `out_point_five$name`, which has elements like "beta[1]+beta[2]",
    ## to fitting windows (rows of `fr`). If `link = FALSE`, then we
    ## also need to inverse link transform `estimate`, `lower`, `upper`.

    fnl <- object$tmb_args$data$fnl
    ffr <- object$tmb_args$data$ffr[!duplicated(object$index), , drop = FALSE]
    d0 <- data.frame(
      name = unlist(lapply(split(names(object$par), rep.int(seq_along(fnl), fnl)), decontrast_beta_names)),
      response = rep.int(pn, fnl),
      term = rep.int(names(fnl), fnl),
      level = unlist(lapply(ffr[names(fnl)], levels)),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    dl <- split(out_point_five, d0$response[match(out_point_five$name, d0$name)])

    g <- function(d, s) {
      ## FIXME: All a bit of a circus...
      i_d0  <- match(d$name, d0$name)
      j_ffr <- match(d0$term[i_d0[1L]], names(ffr))
      i_d   <- match(ffr[[j_ffr]], d0$level[i_d0])
      if (!link) {
        ## FIXME: Assuming here that inverse link function
        ## is increasing and plays nicely with data frames
        d[2:4] <- get_inverse_link(extract_link_string(s))(d[2:4])
        s <- remove_link_string(s)
      }
      cbind(
        par = factor(s),
        fr,
        d[i_d, -1L, drop = FALSE]
      )
    }
    out <- do.call(rbind, Map(g, d = dl[parm], s = parm))
  }

  j_elu <- length(out) - 2:0 # index of "estimate", "lower", "upper"
  if (any(c("R0", "tdoubling") %in% parm0)) {
    s <- if (link) "log_r" else "r"
    delink <- if (link) exp else identity

    if ("R0" %in% parm0) {
      d <- out[out$par == s, , drop = FALSE]
      d[j_elu] <- lapply(delink(d[j_elu]), compute_R0,
                         breaks = breaks, probs = probs)
      d$par <- factor("R0")
      out <- rbind(out, d)
    }
    if ("tdoubling" %in% parm0) {
      d <- out[out$par == s, , drop = FALSE]
      j_eul <- j_elu[c(1L, 3L, 2L)]
      d[j_elu] <- log(2) / delink(d[j_eul])
      d$par <- factor("tdoubling")
      out <- rbind(out, d)
    }

    out$par <- factor(out$par, levels = if (link) parm0 else remove_link_string(parm0))
    out <- out[!is.na(out$par), ]
    out <- out[order(out$par), ]
  }
  row.names(out) <- NULL


  ## For plot.egf_confint()
  refdate <- min(object$frame[[1L]], attr(object$frame, "extra")[[1L]])
  time_split <- split(days(object$frame[[1L]], since = refdate), object$index)

  structure(out,
    level = level,
    refdate = refdate,
    endpoints = data.frame(
      .t1 = vapply(time_split, min, 0L),
      .t2 = vapply(time_split, max, 0L)
    ),
    group_by = attr(object$frame, "group_by"),
    class = c("egf_confint", "data.frame")
  )
}

#' Plot confidence intervals
#'
#' A method for graphically comparing confidence intervals
#' on fitted values across groups.
#'
#' @param x
#'   An `"egf_confint"` object returned by [confint.egf()].
#' @param type
#'   A character string indicating how confidence intervals are
#'   displayed (see Details).
#' @param group_by
#'   A formula of the form `~f1:...:fn` specifying an interaction of the
#'   factors in `x` other than `par`. Confidence intervals belonging to
#'   the same level of the interaction are plotted together in the order
#'   set out by `sort`. Use `~1` (the default) to indicate no grouping
#'   or if `x` has no factors other than `par`. Ignored if `type = "2"`.
#' @param sort
#'   A character string indicating how confidence intervals are sorted
#'   within the groups set out by `group_by`. Ignored if `type = "2"`.
#' @param subset
#'   A named list of atomic vectors with elements specifying levels of
#'   factors in `x`. Only the relevant subset of confidence intervals is
#'   plotted. Use `NULL` (the default) to plot all confidence intervals.
#' @param per_plot
#'   A positive integer. One plot will display at most this many
#'   confidence intervals (`type = "1"`) or panels (`type = "2"`).
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' Confidence intervals are displayed in one- or two-dimensional plots,
#' depending on `type`.
#'
#' If `type = "1"`, then the horizontal axis measures the parameter.
#' Confidence intervals are plotted as stacked horizontal line segments.
#'
#' If `type = "2"`, then the vertical axis measures the parameter,
#' while the horizontal axis measures time. Confidence intervals
#' are plotted as `x` by `y` boxes, where `x` is confidence interval
#' length and `y` is fitting window length. In this case, there is
#' one panel per time series.
#'
#' If endpoints of confidence intervals are `NA`, then dashed lines
#' are drawn from the point estimate to the boundary of the plotting
#' region to indicate missingness.
#'
#' @return
#' `NULL` (invisibly).
#'
#' @export
#' @import graphics
plot.egf_confint <- function(x,
                             type = c("1", "2"),
                             group_by = ~1,
                             sort = c("none", "increasing", "decreasing"),
                             subset = NULL,
                             per_plot = 12L, ...) {
  ### Argument validation #################################

  i_factor <- seq_along(x)[-c(1L, length(x) - 2:0)]
  any_factors <- (length(i_factor) > 0L)

  a <- attributes(x)
  type <- match.arg(type)
  if (type == "2") {
    x <- cbind(x, a$endpoints)
  }

  if (any_factors) {
    factor_names <- names(x)[i_factor]
    if (!is.null(subset)) {
      stop_if_not(
        is.list(subset),
        length(subset) > 0L,
        !is.null(names(subset)),
        m = "`subset` must be a named list or NULL."
      )
      stop_if_not(
        vapply(subset, is.atomic, FALSE),
        lengths(subset) > 0L,
        names(subset) %in% c("par", factor_names),
        !duplicated(names(subset)),
        unlist(Map(`%in%`, subset, lapply(x[names(subset)], levels))),
        m = paste0(
          "`subset` must specify levels of factors in\n",
          "`x[-(length(x)-2:0)]`."
        )
      )
      w <- Reduce(`&`, Map(`%in%`, x[names(subset)], subset))
      stop_if_not(
        any(w),
        m = "`subset` does not match any rows of `x`."
      )
      x <- droplevels(x[w, , drop = FALSE])
    }
  }
  stop_if_not_positive_integer(per_plot)

  ### Split, order confidence intervals ###################

  x_split <- split(x, x$par)

  if (type == "1") {
    sort <- match.arg(sort)
    if (any_factors) {
      stop_if_not(
        inherits(group_by, "formula"),
        length(group_by) == 2L,
        grepl("^(1|([[:alnum:]._]+(:[[:alnum:]._]+)*))$", deparse(group_by[[2L]])),
        all.vars(group_by) %in% factor_names,
        m = paste0(
          "`group_by` must be `~1` or a formula of\n",
          "the form `~f1:...:fn`, with `f1`,...,`fn`\n",
          "naming factors in `x[-c(1, length(x)-2:0)]`."
        )
      )
      group_by <- all.vars(group_by)
      factor_names <- c(group_by, setdiff(factor_names, group_by))
    } else {
      group_by <- character(0L)
    }
    x_split <- lapply(x_split, function(d) {
      ord <- switch(sort,
        none = seq_len(nrow(d)),
        increasing = d$estimate,
        decreasing = -d$estimate
      )
      d[do.call(order, c(d[group_by], list(ord))), ]
    })
  }

  if (type == "1") {
    ### Loop over response variables ######################

    op <- par(
      mar = c(3.5, 5, 1.5, 1),
      cex.axis = 0.8,
      cex.lab = 0.9,
      cex.main = 0.9
    )
    on.exit(par(op))

    for (i in seq_along(x_split)) {
      d <- x_split[[i]]
      xlim <- range(d[c("estimate", "lower", "upper")], na.rm = TRUE)
      is_na_lower <- is.na(d$lower)
      is_na_upper <- is.na(d$upper)
      if (any_factors) {
        yax_labels <- as.character(interaction(d[factor_names], drop = TRUE, sep = ":"))
      }

      ### Loop over plots #################################

      j <- 0L
      while (j < nrow(d)) {
        k <- j + seq_len(min(per_plot, nrow(d) - j))
        plot.new()
        plot.window(xlim = xlim, ylim = c(per_plot + 1, 0), yaxs = "i")
        abline(v = axTicks(side = 1), lty = 3, col = "grey75")
        segments(
          x0 = ifelse(is_na_lower[k], par("usr")[1L], d$lower[k]),
          x1 = d$estimate[k],
          y0 = seq_along(k),
          y1 = seq_along(k),
          lty = c(1, 2)[1L + is_na_lower[k]],
          lwd = c(2, 1)[1L + is_na_lower[k]]
        )
        segments(
          x0 = d$estimate[k],
          x1 = ifelse(is_na_upper[k], par("usr")[2L], d$upper[k]),
          y0 = seq_along(k),
          y1 = seq_along(k),
          lty = c(1, 2)[1L + is_na_upper[k]],
          lwd = c(2, 1)[1L + is_na_upper[k]]
        )
        points(
          x = d$estimate[k],
          y = seq_along(k),
          pch = 21,
          bg = "grey75"
        )
        box()
        axis(side = 1,
          tcl = -0.4,
          mgp = c(3, 0.3, 0)
        )
        axis(side = 2,
          at = seq_along(k),
          labels = yax_labels[k],
          tick = FALSE,
          las = 1,
          mgp = c(3, 0.2, 0)
        )
        title(xlab = names(x_split)[i], line = 2)
        main <- sprintf("%.3g%% CI", 100 * a$level)
        if (any_factors) {
          main <- sprintf("%s by %s", main, paste(factor_names, collapse = ":"))
        }
        title(main, line = 0.25, adj = 0)
        j <- j + per_plot
      }
    }
  } else {
    ### Loop over response variables ######################

    op <- par(
      mfrow = c(per_plot, 1),
      mar = c(0, 5, 1.5, 1),
      oma = c(3.5, 0, 0, 0),
      cex.axis = 0.8,
      cex.lab = 0.9,
      cex.main = 0.9
    )
    on.exit(par(op))

    xlim <- range(x_split[[1L]][c(".t1", ".t2")])
    xax_at <- daxis(
      left = xlim[1L],
      right = xlim[2L],
      refdate = a$refdate,
      plot = FALSE
    )

    for (i in seq_along(x_split)) {
      d <- x_split[[i]]
      d_split <- split(d, interaction(d[a$group_by], drop = TRUE, sep = ":"))
      ylim <- range(d[c("estimate", "lower", "upper")], na.rm = TRUE)

      ### Loop over plots #################################

      j <- 0L
      while (j < length(d_split)) {

        ### Loop over panels ##############################

        for (k in j + seq_len(min(per_plot, length(d_split) - j))) {
          dd <- d_split[[k]]
          plot.new()
          plot.window(xlim = xlim, ylim = ylim, xaxs = "i")
          abline(v = xax_at, lty = 3, col = "grey75")

          ### Loop over confidence intervals ##############

          for (l in seq_len(nrow(dd))) {
            t12 <- unlist(dd[l, c(".t1", ".t2")], use.names = FALSE)
            elu <- unlist(dd[l, c("lower", "upper")], use.names = FALSE)
            if ((lna <- is.na(elu[2L]))) {
              elu[2L] <- par("usr")[3L]
            }
            if ((una <- is.na(elu[3L]))) {
              elu[3L] <- par("usr")[4L]
            }
            polygon(
              x = t12[c(1L, 2L, 2L, 1L)],
              y = elu[c(1L, 1L, 2L, 2L)],
              col = if (lna) NA else "grey75",
              border = "grey50",
              lty = if (lna) 2 else 1
            )
            polygon(
              x = t12[c(1L, 2L, 2L, 1L)],
              y = elu[c(1L, 1L, 3L, 3L)],
              col = if (una) NA else "grey75",
              border = "grey50",
              lty = if (una) 2 else 1
            )
            segments(
              x0 = t12[1L],
              y0 = elu[1L],
              x1 = t12[2L],
              y1 = elu[1L],
              col = "grey50",
              lwd = 2,
            )
          } # loop over confidence intervals

          box()
          axis(side = 2, las = 1, mgp = c(3, 0.7, 0))
          title(ylab = names(x_split)[i], line = 2.5)
          title(main = sprintf("%.3g%% CI for %s", 100 * a$level, names(d_split)[k]), line = 0.25, adj = 0)
        } # loop over panels

        daxis(
          left = xlim[1L],
          right = xlim[2L],
          refdate = a$refdate,
          mgp2 = c(0.25, 1.25)
        )
        title(xlab = "date", line = 2, xpd = FALSE)
        j <- j + per_plot
      } # loop over plots
    } # loop over response variables
  }

  invisible(NULL)
}
