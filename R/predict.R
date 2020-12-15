#' @importFrom TMB MakeADFun sdreport
predict.egf <- function(object, time = NULL, se = FALSE, ...) {
  dots <- list(...)
  if (length(dots) > 0L) {
    stop_if_not(
      vapply(dots, is.atomic, logical(1L)),
      lengths(dots) == 1L,
      names(dots) %in% names(object$frame)[-(1:2)],
      !duplicated(names(dots)),
      mapply("%in%", dots, lapply(object$frame[names(dots)], levels)),
      m = paste0(
        "`list(...)` must specify levels of factors in `object$frame`,\n",
        "e.g., f1 = l1, f2 = l2, and so on (one level per factor)."
      )
    )
  }
  if (length(object$frame) == 2L) {
    i <- 1L
  } else {
    d <- object$frame[-(1:2)]
    il <- !duplicated(object$index)
    if (length(dots) > 0L) {
      il <- il & Reduce("&", lapply(names(dots), function(s) d[[s]] == dots[[s]]))
    }
    nts <- sum(il)
    if (nts == 0L) {
      stop("`list(...)` must specify a nonempty interaction\n",
           "of the factors in `object$frame`.")
    } else if (nts > 1L) {
      if (!all(duplicated(d[il, ])[-1L])) {
        stop("`list(...)` must specify a unique interaction\n",
             "of the factors in `object$frame`.")
      }
    }
    i <- which(il)
  }
  stop_if_not(
    is.numeric(time),
    length(time) > 0L,
    m = "`time` must be numeric and have nonzero length."
  )
  stop_if_not(
    !anyNA(time),
    m = "`time` must not have missing values."
  )
  stop_if_not(
    diff(time) > 0,
    m = "`time` must be increasing."
  )
  tr <- range(object$madf_args$data$t[object$index == object$index[i[1L]]])
  warn_if_not(
    time >= tr[1L] & time <= tr[2L],
    m = paste0(
      "There are elements of `time` outside of\n",
      "the fitting window specified by `list(...)`."
    )
  )
  stop_if_not(
    is.logical(se),
    length(se) == 1L,
    !is.na(se),
    m = "`se` must be TRUE or FALSE."
  )

  object$madf_args$data <- within(object$madf_args$data, {
    predict_flag <- 1L
    se_flag <- 1L * se
    t_new <- time
    X_new <- X[i[1L], , drop = FALSE]
    Z_new <- Z[i[1L], , drop = FALSE]
  })
  madf_out_new <- do.call(MakeADFun, object$madf_args)

  if (se) {
    madf_out_new$fn(object$par[object$inr])
    sdr <- summary(sdreport(madf_out_new), select = "report")
    out <- split(as.data.frame(sdr), factor(rownames(sdr)))
    names(out) <- sub("_new", "", names(out))
    out$log_int_inc <- rbind(NA_real_, out$log_int_inc)
    out <- lapply(out, function(x) {
      names(x) <- c("estimate", "se")
      row.names(x) <- NULL
      cbind(time, x)
    })
    class(out) <- c("egf_predict", "list")
  } else {
    madf_report_new <- madf_out_new$report(object$par)
    out <- data.frame(
      time = time,
      log_cum_inc = madf_report_new$log_cum_inc_new,
      log_int_inc = c(NA_real_, madf_report_new$log_int_inc_new)
    )
  }

  attr(out, "refdate") <- object$frame$date[i]
  out
}
