#' @export
coef.egf <- function(object, link = TRUE, ...) {
  stop_if_not_tf(link)
  dots <- list(...)
  if (length(dots) > 0L) {
    stop_if_not(
      vapply(dots, is.atomic, logical(1L)),
      lengths(dots) > 0L,
      names(dots) %in% names(object$frame)[-(1:2)],
      !duplicated(names(dots)),
      unlist(Map("%in%", dots, lapply(object$frame[names(dots)], levels))),
      m = paste0(
        "`list(...)` must specify nonempty levels\n",
        "of factors in `object$frame`."
      )
    )
  }

  ## Keep one row of `frame` (factors only) per fitting window.
  ## No loss of information here, since factors have one level
  ## within fitting windows.
  fr <- object$frame[!duplicated(object$index), -(1:2), drop = FALSE]

  ## Here `Y[i, ]` lists fitted responses in the group specified
  ## by `fr[i, ]`. There are `length(pn)` responses, one for each
  ## mixed effects model.
  pn <- get_par_names(object, link = TRUE)
  Y <- matrix(object$report$Y_short_as_vector$estimate, ncol = length(pn))

  if (link) {
    colnames(Y) <- pn
  } else {
    f <- function(x, s) get_inverse_link(s)(x)
    Y <- mapply(f, x = as.data.frame(Y), s = extract_link_string(pn))
    colnames(Y) <- remove_link_string(pn)
  }

  d <- cbind(fr, Y)

  ## Select user-specified fitting windows
  if (length(dots) > 0L) {
    d <- d[Reduce("&", Map("%in%", d[names(dots)], dots)), ]
  }
  row.names(d) <- NULL
  d
}
