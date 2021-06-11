#' Extract fixed effects covariance matrix
#'
#' Retrieves the covariance matrix of fixed effects coefficients
#' and random effects covariance parameters.
#'
#' @param object
#'   An `"egf"` object returned by [egf()].
#' @param full
#'   A logical scalar. If `TRUE`, then covariances of random effects
#'   covariance parameters are retained.
#' @param cor
#'   A logical scalar. If `TRUE`, then a correlation matrix is returned.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' Contrary to the documentation of generic function [stats::vcov()],
#' `vcov.egf()` returns the covariance matrix corresponding to the
#' output of [fixef.egf()], not [coef.egf()]. Currently, [coef.egf()]
#' is an alias for [fitted.egf()] and so returns fitted values of
#' nonlinear and dispersion model parameters, not mixed effects model
#' coefficients.
#'
#' @return
#' A covariance or correlation matrix inheriting from class
#' `"egf_vcov"`. The number of rows is equal to
#' the length of the `"beta"` component of `object$best` plus
#' (if `full = TRUE`) the length of the `"theta"` component.
#'
#' @export
#' @importFrom stats cov2cor
vcov.egf <- function(object, full = FALSE, cor = FALSE, ...) {
  stop_if_not_true_false(full)
  stop_if_not_true_false(cor)
  if (!inherits)

  nb <- names(object$best)
  if (full) {
    k <- object$nonrandom
  } else {
    k <- grep("^beta\\[", nb)
  }
  V <- object$report$cov[k, k]
  if (cor) {
    V <- cov2cor(V)
  }
  dimnames(V) <- rep_len(list(nb[k]), 2L)
  class(V) <- c("egf_vcov", "matrix", "array")
  V
}
