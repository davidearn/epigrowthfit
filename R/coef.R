coef.egf <- function(object, log = FALSE, ...) {
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
        len = c(1L, NA),
        opt = levels(object$frame[[s]]),
        sprintf("`%s` must be a subset of `levels(object$frame[[%s]]).`", s, s)
      )
    }
  }
  check(log,
    what = "logical",
    len = c(1L, NA),
    fails = is.na,
    "`log` must be TRUE or FALSE."
  )

  Y <- object$madf_out$report(object$par)$Y
  pn <- get_par_names(object$curve, object$distr, object$include_baseline)
  if (log) {
    colnames(Y) <- sprintf("log_%s", pn)
  } else {
    Y <- exp(Y)
    colnames(Y) <- pn
  }

  d <- cbind(object$frame[-(1:2)], as.data.frame(Y))
  d <- d[!duplicated(object$index), ]
  if (length(dots) > 0L) {
    d <- d[Reduce("&", lapply(names(dots), function(s) d[[s]] %in% dots[[s]])), ]
  }
  row.names(d) <- NULL
  d
}
