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
    delink <- if (link) exp else function(x) x

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
  out
}

#' Confidence intervals from likelihood profiles
#'
#' Computes confidence intervals on fixed effect coefficients,
#' log standard deviations of random effects, and linear combinations
#' thereof from their (univariate) likelihood profiles.
#'
#' @param object
#'   An `"egf_profile"` object returned by [profile.egf()].
#' @param parm
#'   Unused argument included for generic consistency.
#' @inheritParams confint.egf
#'
#' @details
#' For each likelihood profile (level of `object$name`),
#' [TMB::confint.tmbprofile()] is called to approximate
#' (via linear interpolation) the two solutions of
#' `deviance(x) = qchisq(level, df = 1)`. These are
#' the lower and upper confidence limits of interest.
#'
#' @return
#' A data frame with character variable `name` equal to
#' `levels(object$name)` and numeric variables `estimate`,
#' `lower`, and `upper` supplying estimates and confidence
#' intervals.
#'
#' @export
#' @importFrom stats qchisq confint
confint.egf_profile <- function(object, parm, level = 0.95, ...) {
  stop_if_not(
    is.numeric(level),
    length(level) == 1L,
    level > 0,
    level < 1,
    m = "`level` must be a number in the interval (0,1)."
  )
  stop_if_not(
    qchisq(level, df = 1) < max(object$deviance, na.rm = TRUE),
    m = paste0(
      "Maximum deviance must exceed `qchisq(level, df = 1)`.\n",
      "Reprofile with higher `max_level` or retry with lower `level`."
    )
  )

  object_split <- split(object, object$name)
  f <- function(d) d[which.min(d[["deviance"]]), "value"]
  estimate <- vapply(object_split, f, 0)
  g <- function(d) {
    d <- d[c("value", "deviance")]
    names(d) <- c("parameter", "value")
    class(d) <- c("tmbprofile", "data.frame")
    confint(d)
  }
  ci <- do.call(rbind, lapply(object_split, g))

  out <- data.frame(
    name = names(object_split),
    estimate,
    ci,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  attr(out, "level") <- level
  out
}

#' Confidence bands on log predicted incidence
#'
#' Computes confidence bands on log predicted incidence assuming
#' asymptotic normality of the estimator.
#'
#' @param object
#'   An `"egf_predict"` object returned by [predict.egf()].
#'   Must supply standard errors on log predicted incidence.
#' @param log
#'   A logical scalar. If `FALSE`, then log predicted incidence and
#'   the corresponding confidence bands are inverse log transformed.
#' @inheritParams confint.egf_profile
#'
#' @details
#' Confidence bands on log predicted incidence are computed as
#' `estimate + c(-1, 1) * q * se`, where `q = qnorm(0.5 * (1 + level))`
#' and `estimate` and `se` are log predicted incidence and the
#' approximate (delta method) standard errors.
#'
#' @return
#' If `log = TRUE`, then `object` but with variable `se` in each listed
#' data frame replaced with two variables `lower` and `upper` supplying
#' confidence bands on log predicted incidence.
#'
#' Otherwise, the same list but with variables `estimate`, `lower`,
#' and `upper` inverse log transformed in each listed data frame.
#' In this case, the listed data frames are named `cum_inc` and
#' `int_inc`, rather than `log_cum_inc` and `log_int_inc`.
#'
#' @export
#' @importFrom stats qnorm
confint.egf_predict <- function(object, parm, level = 0.95,
                                log = TRUE, ...) {
  stop_if_not(
    "se" %in% names(object$log_cum_inc),
    "se" %in% names(object$log_int_inc),
    m = paste0(
      "`object` must supply standard errors.\n",
      "Repeat `predict()` but with argument `se = TRUE`."
    )
  )
  stop_if_not(
    is.numeric(level),
    length(level) == 1L,
    level > 0,
    level < 1,
    m = "`level` must be a number in the interval (0,1)."
  )

  q <- qnorm(0.5 * (1 + level))
  f <- function(d) {
    elu <- d$estimate + outer(d$se, c(0, -1, 1) * q)
    colnames(elu) <- c("estimate", "lower", "upper")
    data.frame(time = d$time, if (log) elu else exp(elu))
  }
  out <- lapply(object, f)
  if (!log) {
    names(out) <- sub("^log_", "", names(out))
  }
  attr(out, "subset") <- attr(object, "subset")
  attr(out, "refdate") <- attr(object, "refdate")
  attr(out, "level") <- level
  out
}
