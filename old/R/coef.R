#' Extract model parameter estimates
#'
#' @description
#' Methods for extracting estimates of model parameters from objects
#' of class "egf_init" or "egf".
#'
#' @param object
#'   An "egf_init" or "egf" object.
#' @param log
#'   A logical scalar. If `TRUE`, then parameter estimates are
#'   log-transformed.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' The method for class "egf_init" returns
#' `object$theta_init` if `log = FALSE` and
#' `object$log_theta_init` if `log = TRUE`.
#'
#' The method for class "egf" returns
#' `object$theta_fit` if `log = FALSE` and
#' `object$log_theta_fit` if `log = TRUE`.
#'
#' @seealso [egf_init()], [egf()]
#' @name coef.egf
NULL

#' @rdname coef.egf
#' @export
coef.egf_init <- function(object, log = FALSE, ...) {
  check(log,
    what = "logical",
    len = 1,
    no = is.na,
    "`log` must be `TRUE` or `FALSE`."
  )

  if (log) {
    object$log_theta_init
  } else {
    object$theta_init
  }
}

#' @rdname coef.egf
#' @export
coef.egf <- function(object, log = FALSE, ...) {
  check(log,
    what = "logical",
    len = 1,
    no = is.na,
    "`log` must be `TRUE` or `FALSE`."
  )

  if (log) {
    object$log_theta_fit
  } else {
    object$theta_fit
  }
}
