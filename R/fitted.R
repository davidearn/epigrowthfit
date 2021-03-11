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
#' @param subset
#'   A list of logical expressions to be evaluated in the combined
#'   model frame `do.call(cbind, object$frame_glmm)`. Fitted values
#'   are returned only for the indexed fitting windows. Use the
#'   default (`NULL`) to retrieve fitted values for all fitting
#'   windows.
#' @param select
#'   A subset of `get_par_names(object, link = TRUE)` naming
#'   nonlinear model parameters for which fitted values should
#'   be returned. Use the default (`NULL`) to retrieve fitted
#'   values for all nonlinear model parameters.
#' @param link
#'   A logical scalar. If `FALSE`, then fitted values are inverse
#'   link-transformed.
#' @param append_frame
#'   A logical scalar. If `TRUE`, then observed values of covariates
#'   are reported alongside fitted values.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' `coef.egf()` is currently an alias for `fitted.egf()`.
#'
#' @return
#' A data frame inheriting from class `"egf_fitted"`, with variables:
#' \item{`.par`}{
#'   Nonlinear model parameter, from `get_par_names(object, link)`.
#' }
#' \item{`.window`}{
#'   Fitting window, from `levels(object$frame_ts$.window)`.
#' }
#' \item{`.value`}{
#'   Fitted value of parameter `.par` in fitting window `.window`.
#' }
#'
#' If `append_frame = TRUE`, then observed values of covariates
#' from the combined model frame `do.call(cbind, object$frame_glmm)`
#' are also reported. (If a variable occurs more than once there,
#' because it appears in multiple model frames, then only the
#' earliest instance is retained. Except in very unusual cases,
#' all instances of the variable are identical, so no information
#' is lost.)
#'
#' @export
fitted.egf <- function(object, subset = NULL, select = NULL,
                       link = TRUE, append_frame = TRUE, ...) {
  frame <- do.call(cbind, object$frame_glmm)
  frame <- frame[!duplicated(names(frame))]
  wl <- levels(object$frame_ts$.window)
  n <- length(wl)
  pn <- get_par_names(object, link = TRUE)
  p <- length(pn)

  if (is.null(subset)) {
    i <- rep.int(TRUE, n)
  } else {
    e <- substitute(subset)
    l <- eval(e, envir = frame, enclos = parent.frame())
    stop_if_not(
      is.list(l),
      vapply(l, is.vector, FALSE, "logical"),
      lengths(l) == n,
      m = sprintf("`subset` must evaluate to a list\nof logical vectors of length %d.", n)
    )
    lr <- Reduce(`&`, l)
    i <- lr & !is.na(lr)
  }
  if (is.null(select)) {
    j <- rep.int(TRUE, p)
  } else {
    stop_if_not(
      is.vector(select, "character"),
      m = "`select` must be a character vector."
    )
    j <- pn %in% select
  }
  stop_if_not_true_false(link)

  ## `Y[i, ]` lists fitted values in fitting window `i`
  Y <- object$report$Y_as_vector$estimate
  dim(Y) <- c(n, p)

  ## Inverse link transform
  if (link) {
    colnames(Y) <- pn
  } else {
    f <- function(x, s) get_inverse_link(s)(x)
    Y <- mapply(f, x = as.data.frame(Y), s = get_link_string(pn))
    colnames(Y) <- remove_link_string(pn)
  }

  ## Construct long format data frame
  d <- data.frame(
    .par = rep(factor(colnames(Y)[j], levels = colnames(Y)), each = sum(i)),
    .window = factor(wl[i], levels = wl), # to be repeated `sum(j)` times
    .value = as.vector(Y[i, j, drop = FALSE]),
    frame[i, append_frame & j, drop = FALSE], # to be repeated `sum(j)` times
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
