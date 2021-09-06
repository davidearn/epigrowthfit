#' Extract fixed effects covariance matrix
#'
#' Computes the covariance matrix of fixed effects coefficients
#' and (optionally) random effects covariance parameters.
#'
#' @param object
#'   An \code{"\link{egf}"} object.
#' @param full
#'   A \link{logical} flag. If \code{TRUE}, then covariances
#'   of random effects covariance parameters are retained.
#' @param cor
#'   A \link{logical} flag. If \code{TRUE}, then a correlation
#'   matrix is returned.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' If \code{object} was constructed by \code{\link{egf}(se = TRUE)},
#' then the full covariance matrix has already been computed and
#' is stored in the \link{list} \code{object$sdreport}.
#' \code{vcov} reuses this matrix to avoid needless recomputation.
#'
#' @section Note:
#' Contrary to the documentation of generic function \code{\link{vcov}},
#' \code{vcov.egf} returns the covariance matrix corresponding
#' to the output of \code{\link{fixef.egf}}, not \code{\link{coef.egf}}.
#' Currently, \code{\link{coef.egf}} is an alias for \code{\link{fitted.egf}}:
#' it returns fitted values of top level nonlinear model parameters,
#' not fixed effects coefficients.
#'
#' @return
#' A covariance or correlation \link{matrix} inheriting from \link{class}
#' \code{"egf_vcov"}. The number of rows is equal to the length of segment
#' \code{beta} of the full parameter vector \code{c(beta, theta, b)} plus
#' (if \code{full = TRUE}) the length of segment \code{theta}.
#'
#' @examples
#' example("egf", "epigrowthfit")
#' vcov(object, full = FALSE, cor = FALSE)
#' vcov(object, full = TRUE,  cor = FALSE)
#' vcov(object, full = TRUE,  cor = TRUE)
#'
#' @export
#' @importFrom stats cov2cor
#' @importFrom TMB sdreport
vcov.egf <- function(object, full = FALSE, cor = FALSE, ...) {
  stop_if_not_true_false(full)
  stop_if_not_true_false(cor)

  if (is.null(object$sdreport)) {
    if (egf_has_random(object)) {
      warning(wrap(
        "Computing a Hessian matrix for a model with random effects, ",
        "which might take a while. To avoid needless recomputation, ",
        "do 'object$sdreport <- try(TMB::sdreport(object$tmb_out))' ",
        "before trying 'vcov'."
      ))
    }
    object$sdreport <- try(TMB::sdreport(object$tmb_out))
  }
  if (inherits(object$sdreport, "try-error")) {
    stop(wrap(
      "Unable to proceed because 'TMB::sdreport(object$tmb_out)' ",
      "throws the following error:\n\n",
      conditionMessage(attr(object$sdreport, "condition")), "\n\n",
      "Retry 'vcov' after diagnosing and refitting."
    ))
  }

  V <- object$sdreport$cov.fixed
  if (!full) {
    k <- (rownames(V) == "beta")
    V <- V[k, k, drop = FALSE]
  }
  if (cor) {
    V <- cov2cor(V)
  }
  dimnames(V)[1:2] <- list(enum_dupl_string(rownames(V)))
  class(V) <- c("egf_vcov", "matrix", "array")
  V
}
