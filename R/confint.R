#' @export
#' @importFrom stats qchisq confint
#' @importFrom parallel mcmapply makePSOCKcluster clusterExport clusterMap stopCluster
confint.egf <- function(object, parm = get_par_names(object), level = 0.95,
                        link = TRUE,
                        method = c("wald", "profile", "uniroot"),
                        grid_len = 12, max_width = 7,
                        parallel = c("serial", "multicore", "snow"),
                        cores = getOption("egf.cores", 2L),
                        breaks = NULL, probs = NULL, ...) {
  ## FIXME: Behaves as expected but isn't pretty...

  curve_is_exp_in_limit <- (object$curve %in% c("exponential", "logistic", "richards"))
  pn <- get_par_names(object, link = TRUE)
  pn0 <- if (curve_is_exp_in_limit) c("R0", "tdoubling", pn) else pn
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
    stop("`parm = \"R0\"` requires non-NULL `breaks` and `probs`.\n",
         "See ?compute_R0.")
  }
  parm0 <- parm <- add_link_string(unique(remove_link_string(parm)))
  parm[parm %in% c("R0", "tdoubling")] <- "log_r"
  parm <- unique(parm)

  stop_if_not(
    is.numeric(level),
    length(level) == 1L,
    level > 0 && level < 1,
    m = "`level` must be a number in the interval (0,1)."
  )
  stop_if_not_tf(link)
  method <- match.arg(method)
  parallel <- match.arg(parallel)
  if (parallel != "serial") {
    stop_if_not_positive_integer(cores)
  }

  fr <- object$frame[!duplicated(object$index), -(1:2), drop = FALSE]

  if (method == "wald") {
    estimate <- matrix(object$report$Y_short_as_vector$estimate, ncol = length(pn))
    se <- matrix(object$report$Y_short_as_vector$se, ncol = length(pn))
    q <- qchisq(level, df = 1)
    f <- function(j) {
      elu <- estimate[, j] + outer(se[, j], c(0, -1, 1) * sqrt(q))
      colnames(elu) <- c("estimate", "lower", "upper")
      elu
    }
    ml <- lapply(match(parm, pn), f)
    nrl <- vapply(ml, nrow, integer(1L))
    if (!link) {
      g <- function(x, s) get_inverse_link(s)(x)
      ml <- Map(g, x = ml, s = extract_link_string(parm))
      parm <- remove_link_string(parm)
    }
    out <- cbind(
      par = factor(rep(parm, nrl), levels = parm),
      do.call(rbind, rep(list(fr), length(parm))),
      do.call(rbind, ml)
    )
  } else {
    rid <- object$tmb_args$data$rid
    stop_if_not(
      colSums(rid[, parm, drop = FALSE]) == 0L,
      m = paste0(
        "Cannot compute likelihood profiles w.r.t. parameters\n",
        "modeled with random effects:\n",
        paste(colnames(rid)[colSums(rid) == 0L], collapse = ", ")
      )
    )

    if (method == "profile") {
      pr <- profile(object,
        parm = parm,
        decontrast = TRUE,
        max_level = level + min(0.01, 0.1 * (1 - level)),
        grid_len = grid_len,
        parallel = parallel,
        cores = cores
      )
      d <- confint(pr, level = level)
    } else if (method == "uniroot") {
      stop_if_not(
        is.numeric(max_width),
        length(max_width) > 0L,
        is.finite(max_width),
        max_width > 0,
        m = "`max_width` must be a numeric vector with positive elements."
      )

      lin_comb <- make_lin_comb_for_parm(object, parm) # see R/profile.R
      lcl <- lapply(seq_len(nrow(lin_comb)), function(i) lin_comb[i, ])
      wl <- rep(max_width, length.out = length(parm))
      ytol <- qchisq(level, df = 1) / 2 # y := diff(nll) = deviance / 2

      f <- function(lc, w) {
        TMB::tmbroot(object$tmb_out, target = ytol, lincomb = lc,
                     sd.range = w)
      }

      if (parallel == "multicore") {
        cil <- mcmapply(f, lc = lcl, w = wl, SIMPLIFY = FALSE,
                        mc.cores = cores)
      } else if (parallel == "snow") {
        cl <- makePSOCKcluster(cores)
        on.exit(stopCluster(cl))
        clusterExport(cl, "ytol", envir = environment())
        cil <- clusterMap(cl, f, lc = lcl, w = wl)
      } else { # "serial"
        cil <- Map(f, lc = lcl, w = wl)
      }

      estimate <- as.vector(lin_comb %*% object$par[object$nonrandom])
      ci <- do.call(rbind, cil)
      colnames(ci) <- c("lower", "upper")
      d <- data.frame(name = rownames(lin_comb), estimate, ci, stringsAsFactors = FALSE)
    }

    ## To create output matching `method = "wald"`, need to
    ## map `d$name`, which has elements like "beta[1]+beta[2]",
    ## to rows of `fr`
    d_split <- split(d, object$par_info$response[match(d$name, object$par_info$name)])
    ffr <- object$tmb_args$data$ffr[!duplicated(object$index), , drop = FALSE]

    g <- function(d, s) {
      i_pin <- match(d$name, object$par_info$name)
      j_ffr <- match(object$par_info$term[i_pin[1L]], names(ffr))
      i_d   <- match(ffr[[j_ffr]], object$par_info$level[i_pin])
      if (!link) {
        ## FIXME: trouble if link is not increasing or if link
        ## calls functions that don't play nicely with data frames
        l <- extract_link_string(s)
        d[2:4] <- get_inverse_link(l)(d[2:4])
        s <- remove_link_string(s)
      }
      cbind(
        par = factor(rep(s, nrow(fr))),
        fr,
        d[i_d, -1L, drop = FALSE]
      )
    }
    out <- do.call(rbind, Map(g, d_split[parm], parm))
  }

  j_elu <- length(out) - 2:0 # index of "estimate", "lower", "upper"
  if (any(c("R0", "tdoubling") %in% parm0)) {
    s <- if (link) "log_r" else "r"
    f <- if (link) exp else function(x) x

    if ("R0" %in% parm0) {
      d <- out[out$par == s, , drop = FALSE]
      d[j_elu] <- lapply(f(d[j_elu]), compute_R0, breaks = breaks, probs = probs)
      d$par <- factor("R0")
      out <- rbind(out, d)
    }
    if ("tdoubling" %in% parm0) {
      d <- out[out$par == s, , drop = FALSE]
      j_eul <- j_elu[c(1L, 3L, 2L)]
      d[j_elu] <- log(2) / f(d[j_eul])
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

#' @importFrom stats qchisq confint
confint.egf_profile <- function(object, parm, level = 0.95, ...) {
  stop_if_not(
    is.numeric(level),
    length(level) == 1L,
    level > 0 && level < 1,
    m = "`level` must be a number in the interval (0,1)."
  )
  stop_if_not(
    qchisq(level, df = 1) < max(object$deviance, na.rm = TRUE),
    m = paste0(
      "Maximum deviance must exceed `qchisq(level, df = 1)`.\n",
      "Reprofile with higher `max_level` or retry with lower `level`."
    )
  )

  object_split <- split(object, factor(object$name, levels = unique(object$name)))
  f <- function(d) d[which.min(d[["deviance"]]), "value"]
  estimate <- vapply(object_split, f, numeric(1L))
  g <- function(d) {
    d <- d[c("value", "deviance")]
    names(d) <- c("parameter", "value")
    class(d) <- c("tmbprofile", "data.frame")
    confint(d)
  }
  ci <- do.call(rbind, lapply(object_split, g))
  name <- as.character(object$name[!duplicated(object$name)])
  out <- data.frame(name, estimate, ci, stringsAsFactors = FALSE)
  row.names(out) <- NULL
  attr(out, "level") <- level
  out
}

#' @export
#' @importFrom stats qchisq
confint.egf_predict <- function(object, parm, level = 0.95, log = TRUE, ...) {
  stop_if_not(
    is.numeric(level),
    length(level) == 1L,
    !is.na(level),
    level > 0 && level < 1,
    m = "`level` must be a number in the interval (0,1)."
  )

  q <- qchisq(level, df = 1)
  f <- function(x) {
    elu <- x$estimate + outer(x$se, c(0, -1, 1) * sqrt(q))
    colnames(elu) <- c("estimate", "lower", "upper")
    data.frame(time = x$time, if (log) elu else exp(elu))
  }
  out <- lapply(object, f)
  if (!log) {
    names(out) <- sub("^log_", "", names(out))
  }
  attr(out, "refdate") <- attr(object, "refdate")
  attr(out, "level") <- level
  out
}
