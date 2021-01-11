#' Extract covariance matrix of fixed effects coefficients
#'
#' Retrieves the covariance or correlation matrix of the fixed effects
#' coefficients and, optionally, the log standard deviations of random
#' effects.
#'
#' @param object
#'   An `"egf"` object returned by [egf()].
#' @param full
#'   A logical scalar. If `TRUE`, then covariances involving log
#'   standard deviations of random effects are retained.
#' @param cor
#'   A logical scalar. If `TRUE`, then correlations are returned
#'   instead of covariances.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' Contrary to the documentation of generic function [stats::vcov()],
#' `vcov.egf()` does _not_ return the covariance matrix corresponding
#' to the output of [coef.egf()], which returns fitted values,
#' not coefficients, of the generalized linear mixed effects models.
#'
#' @return
#' A matrix of covariances or correlations. The number of rows is
#' equal to the length of the `"beta"` component of `object$par`
#' plus, if `full = TRUE`, the length the `"log_sd_b"` component.
#'
#' @export
#' @importFrom TMB sdreport
#' @importFrom stats cov2cor
vcov.egf <- function(object, full = FALSE, cor = FALSE, ...) {
  stop_if_not_tf(full)
  stop_if_not_tf(cor)

  k <- if (full) object$nonrandom else grep("^beta\\[", names(object$par))
  m <- if (cor) cov2cor(object$report$cov[k, k]) else object$report$cov[k, k]
  dimnames(m) <- rep.int(list(names(object$par)[k]), 2L)
  m
}
