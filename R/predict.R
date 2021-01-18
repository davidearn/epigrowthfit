#' Compute predicted incidence
#'
#' Computes predicted values of log cumulative incidnece
#' and log interval incidence given user-specified time points.
#'
#' @param object
#'   An `"egf"` object returned by [egf()].
#' @param subset
#'   A named list of atomic vectors of length 1 specifying exactly
#'   one level for each factor in `object$frame` (and thus a unique
#'   fitting window). Use the default (`NULL`) if and only if
#'   `object$frame` has no factors.
#' @param time
#'   A numeric vector listing increasing time points in days since
#'   the start of the fitting window specified by `subset`.
#' @param se
#'   A logical scalar. If `TRUE`, then standard errors on predicted
#'   values (log scale) are also reported.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' Elements of `subset` (if non-`NULL`) must have the form
#' `factor_name = level_name`, where `factor_name` is
#' the name of a factor in `object$frame` and `level_name`
#' is an element of `levels(object$frame$factor_name)`.
#'
#' @return
#' A list inheriting from class `"egf_predict"`, containing two data
#' frames, `log_cum_inc` and `log_int_inc`, each with numeric variables
#' `time`, `estimate`, and `se` (if `se = TRUE`).
#'
#' `log_cum_inc$estimate[i]` stores log cumulative incidence at
#' `time[i]`. For `i > 1`, `log_int_inc$estimate[i]` stores log
#' incidence during the interval from `time[i-1]` to `time[i]`.
#'
#' @seealso [confint.egf_predict()]
#' @export
#' @importFrom TMB MakeADFun sdreport
predict.egf <- function(object, subset = NULL, time = NULL,
                        se = FALSE, ...) {
  fr <- object$frame[-(1:2)]
  w <- !duplicated(object$index)
  if (length(fr) > 0L || !is.null(subset)) {
    stop_if_not(
      is.list(subset),
      length(subset) > 0L,
      !is.null(names(subset)),
      m = "`subset` must be a named list or NULL."
    )
    stop_if_not(
      length(subset) == length(fr),
      identical(sort(names(subset)), sort(names(fr))),
      vapply(subset, is.atomic, logical(1L)),
      lengths(subset) == 1L,
      mapply("%in%", subset, lapply(fr[names(subset)], levels)),
      m = paste0(
        "`subset` must specify exactly one level\n",
        "for each factor in `object$frame`."
      )
    )
    w <- Reduce("&", Map("==", fr[names(subset)], subset), w)
    stop_if_not(
      any(w),
      m = "`subset` does not match any fitting windows."
    )
  }
  stop_if_not(
    is.numeric(time),
    length(time) >= 2L,
    m = "`time` must be numeric and have length 2 or greater."
  )
  stop_if_not(
    !anyNA(time),
    m = "`time` must not have missing values."
  )
  stop_if_not(
    diff(time) > 0,
    m = "`time` must be increasing."
  )
  stop_if_not_tf(se)

  i <- which(w)
  tr <- range(object$tmb_args$data$t[object$index %in% object$index[i]])
  warn_if_not(
    time >= tr[1L],
    time <= tr[2L],
    m = paste0(
      "There are elements of `time` outside of\n",
      "the fitting window specified by `subset`."
    )
  )

  ## Create the data objects needed to run prediction code
  ## in C++ template
  sparse_X_flag <- Xs <- Xd <- Z <- NULL # R CMD check
  object$tmb_args$data <- within(object$tmb_args$data, {
    predict_flag <- 1L
    se_flag <- 1L * se
    t_new <- time
    if (sparse_X_flag == 1L) {
      ## sparseMatrix() doesn't recycle like matrix()
      Xs_new <- do.call(rbind, rep.int(list(Xs[i[1L], , drop = FALSE]), length(t_new)))
      Xd_new <- Xd
    } else {
      Xd_new <- matrix(Xd[i[1L], , drop = TRUE], nrow = length(t_new), ncol = ncol(Xd), byrow = TRUE)
      Xs_new <- Xs
    }
    if (has_random(object)) {
      Z_new <- do.call(rbind, rep.int(list(Z[i[1L], , drop = FALSE]), length(t_new)))
    } else {
      Z_new <- Z
    }
  })

  tmb_out_new <- do.call(MakeADFun, object$tmb_args)

  if (se) {
    tmb_out_new$fn(object$par[object$nonrandom])
    ssdr <- split_sdreport(sdreport(tmb_out_new))
    out <- ssdr[c("log_cum_inc_new", "log_int_inc_new")]
    out$log_int_inc_new <- rbind(NA_real_, out$log_int_inc_new)
    out <- lapply(out, function(d) cbind(time, d))
    names(out) <- sub("_new$", "", names(out))
  } else {
    r <- tmb_out_new$report(object$par)
    out <- list(
      log_cum_inc = data.frame(time, estimate = r$log_cum_inc_new),
      log_int_inc = data.frame(time, estimate = c(NA_real_, r$log_int_inc_new))
    )
  }
  attr(out, "subset") <- subset
  attr(out, "refdate") <- object$frame$date[i]
  class(out) <- c("egf_predict", "list")
  out
}

#' Confidence bands on predicted incidence
#'
#' Computes confidence bands on predicted incidence assuming
#' asymptotic normality.
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
#' In this case, the list has names `c("cum_inc", "int_inc")`,
#' rather than `c("log_cum_inc", "log_int_inc")`
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
