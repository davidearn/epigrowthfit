#' @importFrom TMB MakeADFun sdreport
predict.egf <- function(object, time = NULL, se = FALSE, ...) {
  stop_if_not_tf(se)
  dots <- list(...)
  if (length(dots) > 0L) {
    stop_if_not(
      vapply(dots, is.atomic, logical(1L)),
      lengths(dots) == 1L,
      names(dots) %in% names(object$frame)[-(1:2)],
      !duplicated(names(dots)),
      mapply("%in%", dots, lapply(object$frame[names(dots)], levels)),
      m = paste0(
        "`list(...)` must specify valid levels\n",
        "of factors in `object$frame` (one level per factor)."
      )
    )
  }

  ## Make sure user has specified a unique group,
  ## if not a unique fitting window
  d <- object$frame[-(1:2)]
  li <- !duplicated(object$index)
  if (length(d) > 0L) {
    if (length(dots) > 0L) {
      li <- li & Reduce("&", Map("==", d[names(dots)], dots))
    }
    stop_if_not(
      any(li),
      duplicated(d[li, , drop = FALSE])[-1L],
      m = paste0(
        "`list(...)` must specify a unique level of\n",
        "`interaction(object$frame[-(1:2)], drop = TRUE)`."
      )
    )
  }
  i <- which(li)

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
  tr <- range(object$tmb_args$data$t[object$index == object$index[i[1L]]])
  warn_if_not(
    time >= tr[1L] & time <= tr[2L],
    m = paste0(
      "There are elements of `time` outside of\n",
      "the fitting window specified by `list(...)`."
    )
  )

  ## Create data objects needed for prediction code
  ## in C++ template to run
  object$tmb_args$data <- within(object$tmb_args$data, {
    predict_flag <- 1L
    se_flag <- 1L * se
    t_new <- time
    if (sparse_X_flag == 1L) {
      ## FIXME: slow? but sparseMatrix() doesn't recycle like matrix()
      Xs_new <- do.call(rbind, rep(list(Xs[i[1L], , drop = FALSE]), length(t_new)))
      Xd_new <- Xd
    } else {
      Xd_new <- matrix(Xd[i[1L], , drop = TRUE], nrow = length(t_new), ncol = ncol(Xd), byrow = TRUE)
      Xs_new <- Xs
    }
    if (has_random(object)) {
      Z_new <- do.call(rbind, rep(list(Z[i[1L], , drop = FALSE]), length(t_new)))
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
    class(out) <- c("egf_predict", "list")
  } else {
    r <- tmb_out_new$report(object$par)
    out <- list(
      log_cum_inc = data.frame(time, estimate = r$log_cum_inc_new),
      log_int_inc = data.frame(time, estimate = c(NA_real_, r$log_int_inc_new))
    )
  }
  attr(out, "refdate") <- object$frame$date[i]
  out
}
