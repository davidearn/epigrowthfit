#' @export
coef.egf <- function(object, link = TRUE, ...) {
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
  stop_if_not_tf(link)

  fr <- object$frame[!is.na(object$index) & !duplicated(object$index), -(1:2), drop = FALSE]
  pn <- get_par_names(object, link = TRUE)
  Y <- matrix(object$report$Y_short_as_vector$estimate, ncol = length(pn))
  if (link) {
    colnames(Y) <- pn
  } else {
    j <- grep("^log_", pn)
    Y[, j] <- exp(Y[, j])
    j <- grep("^logit_", pn)
    Y[, j] <- 1 / (1 + exp(-Y[, j]))
    colnames(Y) <- remove_link_string(pn)
  }

  d <- cbind(fr, Y)
  if (length(dots) > 0L) {
    f <- function(s) d[[s]] %in% dots[[s]]
    d <- d[Reduce("&", lapply(names(dots), f)), ]
  }
  row.names(d) <- NULL
  d
}
