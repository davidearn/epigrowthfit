#' Extract fitted values
#'
#' Extracts the fitted values of incidence model parameters.
#' The fitted value for a given fitting window is obtained
#' by adding
#' (i) the population fitted value computed from the relevant
#' fixed effect coefficients and
#' (ii) the random effect conditional modes relevant to the
#' window (if any).
#'
#' @param object
#'   An `"egf"` object returned by [egf()].
#' @param subset
#'   A named list of atomic vectors with elements specifying levels
#'   of factors in `object$frame` (and thus fitting windows). Use
#'   the default (`NULL`) to retrieve fitted values for every fitting
#'   window or if `object$frame` has no factors.
#' @param link
#'   A logical scalar. If `FALSE`, then fitted values
#'   are inverse link-transformed.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' Elements of `subset` (if non-`NULL`) must have the form
#' `factor_name = level_names`, where `factor_name` is the
#' name of a factor in `object$frame` and `level_names` is
#' a subset of `levels(object$frame$factor_name)`.
#'
#' Currently, `coef.egf()` is an alias for `fitted.egf()`.
#'
#' @return
#' A data frame inheriting from class `"egf_fitted"`, with one
#' row per fitting window. The first `length(object$frame)-2`
#' variables specify levels for each factor in `object$frame`
#' (and thus fitting windows). The remaining variables list
#' fitted values for each element of `get_par_names(object, link)`.
#'
#' @export
fitted.egf <- function(object, subset = NULL, link = TRUE, ...) {
  fr <- object$frame[!duplicated(object$index), -(1:2), drop = FALSE]
  if (length(fr) > 0L && !is.null(subset)) {
    stop_if_not(
      is.list(subset),
      length(subset) > 0L,
      !is.null(names(subset)),
      m = "`subset` must be a named list or NULL."
    )
    stop_if_not(
      vapply(subset, is.atomic, FALSE),
      lengths(subset) > 0L,
      names(subset) %in% names(fr),
      !duplicated(names(subset)),
      unlist(Map(`%in%`, subset, lapply(fr[names(subset)], levels))),
      m = "`subset` must specify levels of factors in `object$frame`."
    )
    w <- Reduce(`&`, Map(`%in%`, fr[names(subset)], subset))
    stop_if_not(
      any(w),
      m = "`subset` does not match any fitting windows."
    )
  }
  stop_if_not_tf(link)

  ## `Y[i, ]` lists fitted values (link scale) in
  ## fitting window `fr[i, ]`
  pn <- get_par_names(object, link = TRUE)
  Y <- object$report$Y_short_as_vector$estimate
  dim(Y) <- c(nrow(fr), length(pn))

  if (link) {
    colnames(Y) <- pn
  } else {
    f <- function(x, s) get_inverse_link(s)(x)
    Y <- mapply(f, x = as.data.frame(Y), s = extract_link_string(pn))
    colnames(Y) <- remove_link_string(pn)
  }

  ## Select user-specified fitting windows
  d <- cbind(fr, Y)
  if (length(fr) > 0L && !is.null(subset)) {
    d <- d[w, , drop = FALSE]
  }
  row.names(d) <- NULL
  class(d) <- c("egf_fitted", "data.frame")
  d
}

#' @rdname fitted.egf
#' @export
coef.egf <- fitted.egf
