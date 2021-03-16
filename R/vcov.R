#' Extract fixed effects covariance matrix
#'
#' Retrieves the covariance or correlation matrix of fixed effects
#' coefficients and, optionally, log standard deviations of random
#' effects coefficients.
#'
#' @param object
#'   An `"egf"` object returned by [egf()].
#' @param full
#'   A logical scalar. If `TRUE`, then covariances involving log
#'   standard deviations of random effects coefficients are retained.
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
#' nonlinear model parameters, not mixed effects model coefficients.
#'
#' @return
#' A covariance or correlation matrix inheriting from class
#' `"egf_vcov"`. The number of rows is equal to the length
#' of the `"beta"` component of `object$best` plus the length
#' the `"log_sd_b"` component (if `full = TRUE`).
#'
#' @export
#' @importFrom stats cov2cor
vcov.egf <- function(object, full = FALSE, cor = FALSE, ...) {
  stop_if_not_true_false(full)
  stop_if_not_true_false(cor)

  if (full) {
    k <- object$nonrandom
  } else {
    k <- grep("^beta\\[", names(object$best))
  }
  if (cor) {
    m <- cov2cor(object$report$cov[k, k])
  } else {
    m <- object$report$cov[k, k]
  }
  dimnames(m) <- rep_len(list(names(object$best)[k]), 2L)
  class(m) <- c("egf_vcov", "matrix", "array")
  m
}
