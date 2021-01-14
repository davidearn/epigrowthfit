#' Compute predicted incidence time series
#'
#' Computes predicted values of log cumulative and log interval
#' incidence given user-specified time points.
#'
#' @param object
#'   An `"egf"` object returned by [egf()].
#' @param subset
#'   A named list of atomic vectors of length 1 specifying exactly
#'   one level for each factor in `object$frame` (and thus a unique
#'   parametrization of the incidence model). Use the default (`NULL`)
#'   if (and only if) there are no factors in `object$frame`.
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
#' Elements of `subset` must be named and have the form
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
#' `time[i]`. `log_int_inc$estimate[i]` stores log incidence during
#' the interval from `time[i-1]` to `time[i]`.
#'
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

  ## Create data objects needed to run prediction code in C++ template
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
  attr(out, "refdate") <- object$frame$date[i]
  class(out) <- c("egf_predict", "list")
  out
}
