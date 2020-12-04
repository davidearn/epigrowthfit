coef.egf <- function(object, ...) {
  dots <- list(...)
  if (length(dots) > 0L) {
    check(names(dots),
      opt = names(object$frame)[-(1:2)],
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

  Q <- object$madf_out$report(object$par)$Q
  colnames(Q) <- sprintf("log_%s", colnames(object$madf_data$fid))
  d <- data.frame(object$frame[-(1:2)], as.data.frame(Q))
  d <- do.call(rbind, lapply(split(d, object$index), "[", 1L, ))

  if (length(dots) == 0L) {
    return(d)
  }
  d[Reduce("&", lapply(names(dots), function(s) d[[s]] %in% dots[[s]])), ]
}
