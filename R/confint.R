#' @importFrom
confint.egf <- function(object, parm, level = 0.95, log = FALSE,
                        method = c("profile", "uniroot"),
                        breaks = NULL, probs = NULL, ...) {
  pn <- get_par_names(object$curve, object$distr, object$include_baseline)
  opt <- c("doubling_time", pn)
  if (!is.null(breaks) && !is.null(probs)) {
    opt <- c("R0", opt)
  }
  check(parm,
    what = "character",
    len = 1L,
    opt = opt,
    "`parm` must be one of:\n",
    paste(sprintf("\"%s\"", opt), collapse = ", ")
  )
  check(level,
    what = "numeric",
    len = 1L,
    val = c(0, 1),
    "`level` must be a number in the interval [0,1]."
  )
  check(log,
    what = "logical",
    len = 1L,
    fails = is.na,
    "`log` must be TRUE or FALSE."
  )
  method <- match.arg(method)

  if (parm %in% c("doubling_time", "R0")) {
    ci <- confint(object, parm = "r", level = level, method = method, ...)
    i_elu <- length(ci) - 2:0
    ci[i_elu] <- switch(parm,
      doubling_time = lapply(ci[i_elu[c(1L, 3L, 2L)]], compute_doubling_time),
      R0 = lapply(ci[i_elu], compute_R0, breaks, probs)
    )
    return(ci)
  }

  ## Much ado about finding coefficients
  ind <- which(!duplicated(object$index))
  lincomb <- with(object$madf_args$data[c("X", "Z", "fnl", "rnl", "rid")], {
    lc <- matrix(0, nrow = length(ind), ncol = length(object$par))
    k <- match(parm, pn)
    j <- seq(to = cumsum(fnl)[k], length.out = fnl[k])
    lc[, j] <- X[ind, j, drop = FALSE]
    rnp <- rowSums(rid)
    b <- tail(object$par, sum(rnl * rnp))
    v <- rep(seq_along(pn), nrow(rid))
    u <- as.vector(t(rid))
    j <- sum(fnl) + which(v[u > 0L] == k)
    u <- as.vector(t(rnl * rid))
    b <- b[rep(v, u) == k]
    jj <- rep(rid[, k], rnl) > 0L
    lc[, j] <- t(apply(Z[ind, jj, drop = FALSE], 1L, function(x) b[x > 0L]))
    lc
  })

  q <- qchisq(level, df = 1)
  out <- data.frame(
    estimate = as.vector(lincomb %*% object$par),
    lower = NA_real_,
    upper = NA_real_
  )
  if (!log) {
    out$estimate <- exp(out$estimate)
  }

  for (i in seq_along(ind)) {
    if (method == "profile") {
      pf <- tmbprofile(object$madf_out, lincomb = lincomb[i, ], ...)
      log_lu <- as.numeric(confint(pf, level = level))
    } else if (method == "uniroot") {
      log_lu <- unname(tmbroot(object$madf_out, target = 0.5 * q, lincomb = lincomb[i, ], ...))
      # } else if (method == "wald") {
      #   sdr <- summary(sdreport(object$madf_out), select = "fixed")
      #   ese <- sdr[log_parm, c("Estimate", "Std. Error")]
      #   log_lu <- ese[1] + c(-1, 1) * sqrt(q) * ese[2]
    }
    out[i, c("lower", "upper")] <- if (log) log_lu else exp(log_lu)
  }

  out <- cbind(object$frame[ind, -(1:2), drop = FALSE], out)
  row.names(out) <- NULL
  out
}
