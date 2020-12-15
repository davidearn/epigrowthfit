coef.egf <- function(object, log = TRUE, ...) {
  dots <- list(...)
  if (length(dots) > 0L) {
    stop_if_not(
      vapply(dots, is.atomic, logical(1L)),
      lengths(dots) > 0L,
      names(dots) %in% names(object$frame)[-(1:2)],
      !duplicated(names(dots)),
      unlist(Map("%in%", dots, lapply(object$frame[names(dots)], levels))),
      m = paste0(
        "`list(...)` must specify valid levels\n",
        "of factors in `object$frame`."
      )
    )
  }
  stop_if_not(
    is.logical(log),
    length(log) == 1L,
    !is.na(log),
    m = "`log` must be TRUE or FALSE."
  )

  pn <- get_par_names(object)
  sdr <- summary(object$madf_sdreport, select = "report")
  Y <- matrix(sdr[rownames(sdr) == "y", 1L], ncol = length(pn))
  if (log) {
    colnames(Y) <- pn
  } else {
    Y <- exp(Y)
    colnames(Y) <- sub("^log_", "", pn)
  }
  d <- cbind(
    object$frame[!duplicated(object$index), -(1:2), drop = FALSE],
    as.data.frame(Y)
  )
  if (length(dots) > 0L) {
    f <- function(s) d[[s]] %in% dots[[s]]
    d <- d[Reduce("&", lapply(names(dots), f)), ]
  }
  row.names(d) <- NULL
  d
}
