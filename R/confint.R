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
  stop_if_not(
    is.numeric(level),
    length(level) == 1L,
    level > 0,
    level < 1,
    m = "`level` must be a number in the interval (0,1)."
  )
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
    out <- out[order(out$par), ]
    out <- out[!is.na(out$par), ]
  }

  row.names(out) <- NULL
  attr(out, "level") <- level
  class(out) <- c("egf_confint", "data.frame")
  out
}

#' Plot confidence intervals
#'
#' A method for graphically comparing confidence intervals
#' on fitted values across groups.
#'
#' @param x
#'   An `"egf_confint"` object returned by [confint.egf()].
#' @param par
#'   A subset of `levels(x$par)` specifying response variables
#'   for which confidence intervals are desired.
#' @param group_by
#'   A formula of the form `~f1:...:fn` specifying an interaction of
#'   the factors in `x` other than `par`. `group_by` determines the
#'   order in which confidence intervals are plotted, with `f1` varying
#'   slowest. Use `~1` (the default) to preserve the order given in `x`
#'   or if `x` has no factors other than `par`.
#' @param subset
#'   A named list of atomic vectors with elements specifying
#'   levels of factors in `x` other than `par`. Only the subset
#'   of confidence intervals belonging to these levels is plotted.
#'   Use `NULL` (the default) to plot all confidence intervals.
#' @param ci_per_panel
#'   A positive integer. One panel will display at most this many
#'   confidence intervals.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' `NULL` (invisibly).
#'
#' @export
#' @import graphics
plot.egf_confint <- function(x, par = levels(d$par),
                             group_by = ~1, subset = NULL,
                             ci_per_panel = 8L, ...) {
  stop_if_not(
    is.atomic(par),
    length(par) > 0L,
    par %in% levels(x$par),
    !duplicated(par),
    m = "`par` must be a subset of `levels(x$par)`."
  )
  any_factors <- (length(x) > 4L)
  if (any_factors) {
    factor_names <- names(x)[-c(1L, length(x) - 2:0)]
    stop_if_not(
      inherits(group_by, "formula"),
      length(group_by) == 2L,
      grepl("^(1|([[:alnum:]._]+(:[[:alnum:]._]+)*))$", deparse(group_by[[2L]])),
      all.vars(group_by) %in% factor_names,
      m = paste0(
        "`group_by` must be a formula of the form\n",
        "`~1` or `~f1:...:fn`, with `f1`,...,`fn`\n",
        "naming factors in `x[-c(1, length(x)-2:0)]`."
      )
    )
    group_by <- all.vars(group_by)
    any_groups <- (length(group_by) > 0L)

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
        names(subset) %in% factor_names,
        !duplicated(names(subset)),
        unlist(Map("%in%", subset, lapply(x[names(subset)], levels))),
        m = paste0(
          "`subset` must specify levels of factors in",
          "`x[-c(1, length(x)-2:0)]`."
        )
      )
      w <- Reduce("&", Map("%in%", x[names(subset)], subset))
      stop_if_not(
        any(w),
        m = "`subset` does not match any rows of `x`."
      )
      x <- droplevels(x[w, , drop = FALSE])
      factor_names <- c(group_by, setdiff(factor_names, group_by))
    }
  } else {
    any_groups <- FALSE
  }
  stop_if_not_positive_integer(ci_per_panel)

  x_split <- split(x, factor(x$par, levels = par))
  if (any_groups) {
    x_split <- lapply(x_split, function(d) d[do.call(order, d[group_by]), ])
  }

  op <- graphics::par(
    mar = c(3.5, 5, 1.5, 1),
    cex.axis = 0.8,
    cex.lab = 0.9,
    cex.main = 0.9
  )
  on.exit(graphics::par(op))

  for (i in seq_along(x_split)) {
    d <- x_split[[i]]
    is_na_lower <- is.na(d$lower)
    is_na_upper <- is.na(d$upper)
    if (any_factors) {
      yax_labels <- as.character(interaction(d[factor_names], drop = TRUE, sep = ":"))
    }

    j <- 0L
    while (j < nrow(d)) {
      k <- j + seq_len(min(ci_per_panel, nrow(d) - j))
      plot.new()
      plot.window(
        xlim = range(d[c("estimate", "lower", "upper")], na.rm = TRUE),
        ylim = c(ci_per_panel + 1, 0),
        yaxs = "i"
      )
      abline(v = axTicks(side = 1), lty = 3, col = "grey75")
      segments(
        x0 = ifelse(is_na_lower[k], graphics::par("usr")[1L], d$lower[k]),
        x1 = d$estimate[k],
        y0 = seq_along(k),
        y1 = seq_along(k),
        lty = c(1, 3)[is_na_lower[k] + 1L],
        lwd = c(2, 1)[is_na_lower[k] + 1L]
      )
      segments(
        x0 = d$estimate[k],
        x1 = ifelse(is_na_upper[k], graphics::par("usr")[2L], d$upper[k]),
        y0 = seq_along(k),
        y1 = seq_along(k),
        lty = c(1, 3)[1L + is_na_upper[k]],
        lwd = c(2, 1)[1L + is_na_upper[k]]
      )
      points(
        x = d$estimate[k],
        y = seq_along(k),
        pch = 21,
        bg = "grey75"
      )
      box()
      axis(side = 1, tcl = -0.4, mgp = c(3, 0.3, 0))
      axis(side = 2,
        at = seq_along(k),
        labels = yax_labels[k],
        tick = FALSE,
        las = 1,
        mgp = c(3, 0.2, 0)
      )
      xlab <- names(x_split)[i]
      if (grepl("_", xlab)) {
        ss <- strsplit(xlab, "_")[[1L]]
        xlab <- sprintf("%s(%s)", ss[1L], ss[2L])
      }
      title(xlab = xlab, line = 2)
      main <- sprintf("%.3g%% CI", 100 * attr(x, "level"))
      if (any_factors) {
        main <- sprintf("%s by %s", main, paste(factor_names, collapse = ":"))
      }
      title(main, line = 0.25, adj = 0)
      j <- j + ci_per_panel
    }
  }
  invisible(NULL)
}
