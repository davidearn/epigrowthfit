#' Extract fitted values
#'
#' Extracts the fitted values of nonlinear model parameters.
#' The fitted value for a given fitting window is obtained
#' by adding
#' (i) the population fitted value computed from the relevant
#' fixed effects coefficients and
#' (ii) all random effects terms (if any), with random effects
#' coefficients set equal to their conditional modes.
#'
#' @param object
#'   An `"egf"` object returned by [egf()].
#' @param par
#'   A subset of `get_par_names(object, link = TRUE)` naming nonlinear
#'   model parameters for which fitted values should be retrieved.
#' @param subset
#'   An expression to be evaluated in the combined model frame. Must
#'   evaluate to a logical vector or list of logical vectors indexing
#'   rows of the data frame, and thus fitting windows. Fitted values
#'   are retrieved only for the indexed fitting windows. The default
#'   (`NULL`) is to consider all fitting windows.
#' @param append
#'   An expression indicating variables in the combined model frame
#'   to be included with the result. The default (`NULL`) is to append
#'   nothing.
#' @param link
#'   A logical scalar. If `FALSE`, then fitted values are inverse
#'   link-transformed.
#' @param .subset
#'   A logical vector, to be used (if non-`NULL`) in place
#'   of the result of evaluating `subset`.
#' @param .append
#'   A character vector listing variable names, to be used
#'   (if non-`NULL`) in place of the result of evaluating `append`.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' `coef.egf()` is currently an alias for `fitted.egf()`.
#'
#' # Nonstandard evaluation
#'
#' `subset` and `append` are evaluated in a nonstandard way
#' to make interactive use more convenient. They are handled
#' much like the [subset()] arguments `subset` and `select`.
#' To avoid unexpected behaviour, especially when programming,
#' use `.subset` and `.append`.
#'
#' # Warning
#'
#' The combined model frame is
#' `do.call(cbind, unname(object$frame_par))`.
#' If a variable occurs more than once there, because it appears
#' in multiple model frames, then only the earliest instance is
#' retained. Except in unusual cases, all instances of a variable
#' will be identical, so no information will be lost.
#'
#' @return
#' A data frame inheriting from class `"egf_fitted"`, with variables:
#' \item{`par`}{
#'   Nonlinear model parameter,
#'   from `get_par_names(object, link = link)`.
#' }
#' \item{`ts`}{
#'   Time series, from `levels(object$frame_ts$ts)`.
#' }
#' \item{`window`}{
#'   Fitting window, from `levels(object$frame_ts$window)`.
#' }
#' \item{`value`}{
#'   Fitted value of parameter `par` in fitting window `window`.
#' }
#'
#' @export
fitted.egf <- function(object,
                       par = get_par_names(object, link = TRUE),
                       subset = NULL,
                       append = NULL,
                       link = TRUE,
                       .subset = NULL,
                       .append = NULL,
                       ...) {
  stop_if_not_true_false(link)

  frame <- do.call(cbind, unname(object$frame_par))
  frame <- frame[!duplicated(names(frame))]

  par <- unique(match.arg(par, several.ok = TRUE))
  subset <- subset_to_index(substitute(subset), frame, parent.frame(),
                            .subset = .subset)
  append <- append_to_index(substitute(append), frame, parent.frame(),
                            .append = .append)

  ts <- object$frame_ts$ts
  window <- object$frame_ts$window
  k <- !is.na(window) & !duplicated(window)
  pn <- get_par_names(object, link = TRUE)

  ## `Y[i, j]` is fitted value of nonlinear model parameter `j`
  ## in fitting window `i`
  Y <- object$report$Y_as_vector$estimate
  dim(Y) <- c(nlevels(window), length(pn))
  colnames(Y) <- pn

  Y <- Y[subset, par, drop = FALSE]
  if (link) {
    f <- identity
  } else {
    Y <- mapply(function(x, s) get_inverse_link(s)(x),
      x = as.data.frame(Y),
      s = get_link_string(par)
    )
    f <- remove_link_string
  }

  d <- data.frame(
    par = rep(factor(par, levels = pn, labels = f(pn)), each = length(subset)),
    ts = ts[k][subset],
    window = window[k][subset],
    value = as.numeric(Y),
    frame[subset, append, drop = FALSE],
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  class(d) <- c("egf_fitted", "data.frame")
  d
}

#' @rdname fitted.egf
#' @export
coef.egf <- fitted.egf
