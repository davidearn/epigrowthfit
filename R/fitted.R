#' Extract fitted values
#'
#' Extracts the fitted values of incidence model parameters.
#' The fitted value for a given fitting window is obtained
#' by adding
#' (i) the population fitted value computed from the relevant
#' fixed effect coefficients and
#' (ii) the random effect conditional modes relevant to the
#' window.
#'
#' @param object
#'   An `"egf"` object returned by [egf()].
#' @param link
#'   A logical scalar. If `FALSE`, then inverse link-transformed
#'   fitted values are returned.
#' @param ...
#'   Atomic vectors, optionally specifying fitting windows
#'   for which fitted values are desired. See Details.
#'
#' @details
#' Elements of `list(...)` must be named and have the form
#' `factor_name = c(level_names)`, where `factor_name` must
#' be the name of a factor in `object$frame` and `level_names`
#' must be a subset of `levels(object$frame[[factor_name]])`.
#'
#' @return
#' A data frame inheriting from class `"egf_fitted"`, with one row
#' per fitting window. The first `length(object$frame)-2` variables
#' specify levels for each factor in `object$frame[-(1:2)]`
#' (and thus fitting windows). The remaining variables list fitted
#' values for each element of `get_par_names(object, link)`.
#'
#' @export
fitted.egf <- function(object, link = TRUE, ...) {
  stop_if_not_tf(link)
  dots <- list(...)
  if (length(dots) > 0L) {
    stop_if_not(
      vapply(dots, is.atomic, FALSE),
      lengths(dots) > 0L,
      names(dots) %in% names(object$frame)[-(1:2)],
      !duplicated(names(dots)),
      unlist(Map("%in%", dots, lapply(object$frame[names(dots)], levels))),
      m = "`list(...)` must specify levels of factors in `object$frame`."
    )
  }

  ## Keep one row of `frame` (factors only) per fitting window.
  ## No loss of information here, since factors have one level
  ## within fitting windows.
  fr <- object$frame[!duplicated(object$index), -(1:2), drop = FALSE]

  ## `Y[i, ]` lists fitted responses in the group specified by `fr[i, ]`
  pn <- get_par_names(object, link = TRUE)
  Y <- matrix(object$report$Y_short_as_vector$estimate, ncol = length(pn))

  if (link) {
    colnames(Y) <- pn
  } else {
    f <- function(x, s) get_inverse_link(s)(x)
    Y <- mapply(f, x = as.data.frame(Y), s = extract_link_string(pn))
    colnames(Y) <- remove_link_string(pn)
  }

  ## Select user-specified fitting windows
  d <- cbind(fr, Y)
  if (length(dots) > 0L) {
    d <- d[Reduce("&", Map("%in%", d[names(dots)], dots)), , drop = FALSE]
  }
  row.names(d) <- NULL
  class(d) <- c("egf_fitted", "data.frame")
  d
}

#' @rdname fitted.egf
#' @export
coef.egf <- fitted.egf
