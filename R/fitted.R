#' Extract fitted values
#'
#' Extracts fitted values of nonlinear model parameters.
#' The fitted value for a given fitting window is obtained
#' by adding
#' (i) the population fitted value computed from the
#' relevant fixed effects coefficients and
#' (ii) all (if any) random effects, with random effects
#' coefficients set equal to their conditional modes.
#'
#' @param object
#'   An `"egf"` object returned by [egf()].
#' @param par
#'   A subset of `get_par_names(object, link = TRUE)` naming nonlinear
#'   model parameters for which fitted values should be retrieved.
#' @param subset
#'   An expression to be evaluated in the combined model frame
#'   (see [make_combined()]). Must evaluate to a logical vector
#'   indexing rows of the data frame, and thus fitting windows.
#'   Fitted values are retrieved only for indexed windows.
#'   The default (`NULL`) is to consider all windows.
#' @param append
#'   An expression indicating variables in the combined model frame
#'   (see [make_combined()]) to be included with the result.
#'   The default (`NULL`) is to append nothing.
#' @param link
#'   A logical scalar. If `FALSE`, then fitted values are inverse
#'   link-transformed.
#' @param se
#'   A logical scalar. If `link = TRUE` and `se = TRUE`, then
#'   approximate (delta method) standard errors on fitted values
#'   are reported. Note that standard errors are required for
#'   subsequent use of [confint.egf_fitted()].
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
#' See topic [`nse`] for details on nonstandard evaluation
#' of `subset` and `append`.
#'
#' @return
#' A data frame inheriting from class `"egf_fitted"`, with variables:
#' \item{`par`}{
#'   Nonlinear model parameter,
#'   from `get_par_names(object, link = link)`.
#' }
#' \item{`ts`}{
#'   Time series, from `levels(object$endpoints$ts)`.
#' }
#' \item{`window`}{
#'   Fitting window, from `levels(object$endpoints$window)`.
#' }
#' \item{`estimate`}{
#'   Fitted value of parameter `par` in fitting window `window`.
#' }
#' \item{`se`}{
#'   (If `link = TRUE` and `se = TRUE`.)
#'   Approximate (delta method) standard error on `value`.
#' }
#'
#' @seealso [confint.egf_fitted()]
#' @export
fitted.egf <- function(object,
                       par = get_par_names(object, link = TRUE),
                       subset = NULL,
                       append = NULL,
                       link = TRUE,
                       se = FALSE,
                       .subset = NULL,
                       .append = NULL,
                       ...) {
  stop_if_not_true_false(link)

  combined <- make_combined(object)
  pn <- get_par_names(object, link = TRUE)
  par <- unique(match.arg(par, several.ok = TRUE))
  subset <- subset_to_index(substitute(subset), combined, parent.frame(),
                            .subset = .subset)
  append <- append_to_index(substitute(append), combined, parent.frame(),
                            .append = .append)

  ## `Y[i, j]` is fitted value of nonlinear model parameter `j`
  ## in fitting window `i`. `Y_se[i, j]` is the standard error.
  Y <- object$report$Y_as_vector$estimate
  Y_se <- object$report$Y_as_vector$se
  dim(Y) <- dim(Y_se) <- c(nrow(object$endpoints), length(pn))
  colnames(Y) <- colnames(Y_se) <- pn
  Y <- Y[subset, par, drop = FALSE]
  Y_se <- Y_se[subset, par, drop = FALSE]

  d <- data.frame(
    par = rep(factor(par, levels = pn), each = length(subset)),
    ts = object$endpoints$ts[subset],
    window = object$endpoints$window[subset],
    estimate = as.numeric(Y),
    combined[subset, append, drop = FALSE],
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (link && se) {
    d$se <- as.numeric(Y_se)
  }
  if (!link) {
    d$estimate <- apply_inverse_link(d$estimate, g = d$par)
  }
  attr(d, "se") <- link && se
  class(d) <- c("egf_fitted", "data.frame")
  d
}

#' @rdname fitted.egf
#' @export
coef.egf <- fitted.egf

#' Confidence intervals on fitted values
#'
#' Computes confidence intervals on fitted values of nonlinear
#' model parameters.
#'
#' @param object
#'   An `"egf_fitted"` object returned by [fitted.egf()].
#'   Must supply standard errors on fitted values.
#' @inheritParams confint.egf
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' Confidence limits on fitted values (link scale) are computed
#' as `estimate + c(-1, 1) * sqrt(q) * se`, with `estimate` and
#' `se` obtained from `object` and `q = qchisq(level, df = 1)`.
#'
#' @return
#' If `link = TRUE`, then `object` but with variable `se` replaced
#' with variables `lower` and `upper` supplying confidence limits
#' on link scale fitted values.
#'
#' Otherwise, the same object but with variables `estimate`, `lower`,
#' and `upper` inverse link transformed and link prefixes stripped
#' from `levels(par)`.
#'
#' `level` is retained as an attribute.
#'
#' @export
confint.egf_fitted <- function(object, parm, level = 0.95, link = TRUE, ...) {
  stop_if_not(
    attr(object, "se"),
    m = "`object` must supply standard errors on fitted values.\nRepeat `fitted()` with `link = TRUE` and `se = TRUE`."
  )
  stop_if_not_number_in_interval(level, 0, 1, "()")
  stop_if_not_true_false(link)

  s <- c("par", "ts", "window", "estimate", "se")
  d <- data.frame(
    object[s[1:4]],
    do_wald(estimate = object$estimate, se = object$se, level = level),
    object[-match(s, names(object), 0L)],
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  attr(d, "level") <- level
  if (link) {
    return(d)
  }
  s_elu <- c("estimate", "lower", "upper")
  d[s_elu] <- apply_inverse_link(d[s_elu], g = d$par)
  levels(d$par) <- remove_link_string(levels(d$par))
  d
}
