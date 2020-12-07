#' @importFrom TMB MakeADFun sdreport
predict.egf <- function(object, time = NULL, se = FALSE, ...) {
  dots <- list(...)
  if (length(dots) > 0L) {
    check(dots,
      each_passes = is.atomic,
      "Elements of `list(...)` must be atomic vectors."
    )
    check(names(dots),
      opt = names(object$frame)[-(1:2)],
      fails = function(x) any(duplicated(x)),
      "`names(list(...))` must name factors in `object$frame`."
    )
    for (s in names(dots)) {
      check(dots[[s]],
        len = 1L,
        opt = levels(object$frame[[s]]),
        sprintf("`%s` must be an element of `levels(object$frame[[%s]]).`", s, s)
      )
    }
  }

  d <- object$frame[-(1:2)]
  if (ncol(d) == 0L) {
    i <- 1L
  } else {
    li <- !duplicated(object$index)
    if (length(dots) > 0L) {
      li <- li & Reduce("&", lapply(names(dots), function(s) d[[s]] == dots[[s]]))
    }
    i <- which(li)
    if (length(i) == 0L) {
      stop("`list(...)` must specify a nonempty interaction.")
    }
    if (length(i) > 1L && sum(!duplicated(d[i, ])) > 1L) {
      stop("`list(...)` must specify a unique interaction.")
    }
  }

  check(time,
    what = "numeric",
    len = c(1L, NA),
    "`time` must be numeric and have nonzero length."
  )
  check(time,
    fails = anyNA,
    "`time` must not have missing values."
  )
  check(time,
    passes = function(x) all(diff(x) > 0),
    "`time` must be increasing."
  )
  check(time,
    val = range(object$madf_args$data$t[object$index == object$index[i[1L]]]),
    action = "warn",
    "There are elements of `time` outside of the fitting window."
  )
  check(se,
    what = "logical",
    len = 1L,
    fails = is.na,
    "`se` must be TRUE or FALSE."
  )

  object$madf_args$data <- within(object$madf_args$data, {
    predict_flag <- 1L
    se_flag <- 1L * se
    t_new <- time
    X_new <- X[i[1L], , drop = FALSE]
    Z_new <- Z[i[1L], , drop = FALSE]
  })
  madf_out <- do.call(MakeADFun, object$madf_args)

  if (se) {
    madf_out$fn(object$par) # ?!
    sdr <- summary(sdreport(madf_out), select = "report")
    sdr_split <- split(as.data.frame(sdr), factor(rownames(sdr)))
    names(sdr_split) <- sub("_new", "", names(sdr_split))
    sdr_split$log_int_inc <- rbind(NA_real_, sdr_split$log_int_inc)
    out <- lapply(sdr_split, function(x) {
      names(x) <- c("estimate", "se")
      row.names(x) <- NULL
      cbind(time, x)
    })
    class(out) <- c("egf_pred", "list")
  } else {
    r <- madf_out$report(object$par)
    out <- data.frame(
      time = time,
      log_cum_inc = r$log_cum_inc_new,
      log_int_inc = c(NA_real_, r$log_int_inc_new)
    )
  }

  structure(out, refdate = object$madf_args$data$date[i])
}
