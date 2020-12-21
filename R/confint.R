#' @importFrom TMB tmbroot
#' @importFrom stats qchisq confint
#' @importFrom parallel mcmapply makePSOCKcluster clusterMap stopCluster
confint.egf <- function(object, parm = "r", level = 0.95,
                        method = c("wald", "profile", "uniroot"),
                        grid_len = 12, max_width = 7,
                        parallel = c("serial", "multicore", "snow"),
                        cores = getOption("egf.cores", 2L),
                        breaks = NULL, probs = NULL, ...) {
  ## FIXME: Simplify string handling, e.g., "log_r" versus "r".
  pn <- colnames(object$madf_args$data$rid)
  pn0 <- c("R0", "doubling_time", sub("^log_", "", pn))
  stop_if_not(
    is.character(parm),
    length(parm) > 0L,
    parm %in% pn0,
    m = paste0(
      "`parm` must be a subset of:\n",
      paste(sprintf("\"%s\"", pn0), collapse = ", ")
    )
  )
  if ("R0" %in% parm && (is.null(breaks) || is.null(probs))) {
    stop("`parm = \"R0\"` requires non-NULL `breaks` and `probs`.\n",
         "See ?compute_R0.")
  }
  parm0 <- parm
  parm[parm %in% c("R0", "doubling_time")] <- "r"
  parm <- sprintf("log_%s", unique(parm))

  stop_if_not(
    is.numeric(level),
    length(level) == 1L,
    level > 0 && level < 1,
    m = "`level` must be a number in the interval (0,1)."
  )
  method <- match.arg(method)
  parallel <- match.arg(parallel)
  if (parallel != "serial") {
    stop_if_not(
      is.numeric(cores),
      length(cores) == 1L,
      is.finite(cores),
      cores >= 1,
      m = "`cores` must be 1 or greater."
    )
  }

  if (method == "wald") {
    fr <- object$frame[!duplicated(object$index), -(1:2), drop = FALSE]
    estimate <- matrix(object$report$y$estimate, ncol = length(pn))
    se <- matrix(object$report$y$se, ncol = length(pn))
    q <- qchisq(level, df = 1)
    f <- function(j) {
      elu <- estimate[, j] + outer(se[, j], c(0, -1, 1) * sqrt(q))
      colnames(elu) <- c("estimate", "lower", "upper")
      elu
    }
    ml <- lapply(match(parm, pn), f)
    out <- cbind(
      par = factor(rep(parm, each = nrow(fr))),
      do.call(rbind, rep(list(fr), length(parm))),
      do.call(rbind, ml)
    )
  } else {
    rid <- object$madf_args$data$rid
    stop_if_not(
      colSums(rid[, parm, drop = FALSE]) == 0L,
      m = paste0(
        "Cannot compute likelihood profile w.r.t. parameters\n",
        "modeled with random effects:\n",
        paste(colnames(rid)[colSums(rid) > 0L], collapse = ", "), "\n",
        "Use `method = \"wald\"` instead."
      )
    )

    if (method == "profile") {
      pr <- profile(object,
        parm = sub("^log_", "", parm),
        decontrast = TRUE,
        max_level = level + min(0.01, 0.1 * (1 - level)),
        grid_len = grid_len,
        parallel = parallel,
        cores = cores
      )
      out <- confint(pr, level = level)[-1L]
    } else if (method == "uniroot") {
      stop_if_not(
        is.numeric(max_width),
        is.finite(max_width),
        max_width > 0,
        m = paste0(
          "`max_width` must be a nonempty numeric vector\n",
          "with positive elements."
        )
      )

      lin_comb <- make_lin_comb(object, parm, decontrast = TRUE)
      a <- attributes(lin_comb)
      pin <- droplevels(a$par_info[a$index, ])
      estimate <- as.vector(lin_comb %*% object$par[object$nonrandom])

      target <- 0.5 * qchisq(level, df = 1)
      lcl <- lapply(seq_len(nrow(lin_comb)), function(i) lin_comb[i, ])
      wl <- rep(max_width, length.out = length(parm))
      f <- function(lc, w) {
        tmbroot(object$madf_out, target = target, lincomb = lc, sd.range = w)
      }
      if (parallel == "multicore") {
        ci <- mcmapply(f, lc = lcl, w = wl, SIMPLIFY = FALSE,
                       mc.cores = cores)
      } else if (parallel == "snow") {
        cl <- makePSOCKcluster(cores)
        ci <- clusterMap(cl, f, lc = lcl, w = wl)
        stopCluster(cl)
      } else { # "serial"
        ci <- Map(f, lc = lcl, w = wl)
      }
      ci <- do.call(rbind, ci)
      colnames(ci) <- c("lower", "upper")
      out <- data.frame(pin, estimate, ci)
    }
  }

  i_elu <- length(out) - 2:0 # index of "estimate", "lower", "upper"
  out[i_elu] <- exp(out[i_elu])
  levels(out$par) <- sub("^log_", "", levels(out$par))

  if ("R0" %in% parm0) {
    d <- out[out$par == "r", ]
    d[i_elu] <- lapply(d[i_elu], compute_R0, breaks = breaks, probs = probs)
    levels(d$par) <- sub("^r$", "R0", levels(d$par))
    out <- rbind(out, d)
  }
  if ("doubling_time" %in% parm0) {
    i_eul <- i_elu[c(1L, 3L, 2L)]
    d <- out[out$par == "r", ]
    d[i_elu] <- log(2) / d[i_eul]
    levels(d$par) <- sub("^r$", "doubling_time", levels(d$par))
    out <- rbind(out, d)
  }
  if ("log_r" %in% parm && !"r" %in% parm0) {
    out <- droplevels(out[out$par != "r", ])
  }
  out <- do.call(rbind, split(out, out$par)[parm0])

  row.names(out) <- NULL
  attr(out, "level") <- level
  out
}

#' @importFrom stats qchisq confint
confint.egf_profile <- function(object, parm, level = 0.95, ...) {
  stop_if_not(
    is.numeric(level),
    length(level) == 1L,
    !is.na(level),
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

  object_split <- split(object, object$name)

  f <- function(d) d[which.min(d[["deviance"]]), "value"]
  estimate <- vapply(object_split, f, numeric(1L))

  g <- function(d) {
    d <- d[c("value", "deviance")]
    names(d) <- c("parameter", "value")
    class(d) <- c("tmbprofile", "data.frame")
    confint(d)
  }
  ci <- do.call(rbind, lapply(object_split, g))

  has_par_info <- ("par" %in% names(object))
  if (has_par_info) {
    pin <- object[!duplicated(object$name), -(length(object) - 1:0), drop = FALSE]
    out <- data.frame(pin, estimate, ci)
    row.names(out) <- NULL
  } else {
    out <- data.frame(estimate, ci)
  }

  attr(out, "level") <- level
  out
}

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
