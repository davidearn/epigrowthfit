#' Update an object returned by an S3 method
#'
#' When \code{\link{match.call}} is called from an S3 method for a generic
#' function, a \link{call} to the method is returned, rather than a call
#' to the generic function. Typically, the method is not exported, so the
#' default method for \code{\link{update}} is unable to evaluate the
#' call in the global environment. These methods for \code{\link{getCall}}
#' circumvent this issue for objects created by various methods implemented
#' in \pkg{epigrowthfit}.
#'
#' @param x
#'   A \link{list} with an element \code{call}, which should be a \link{call}.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' \code{x$call}, modified so that \code{x$call[[1L]]} is a generic function,
#' not an S3 method.
#'
#' @name getCall.egf
NULL

#' @rdname getCall.egf
#' @export
#' @importFrom stats getCall
getCall.egf <- function(x, ...) {
  call <- NextMethod("getCall")
  call[[1L]] <- quote(egf)
  call
}

#' @rdname getCall.egf
#' @export
#' @importFrom stats getCall
getCall.egf_no_fit <- getCall.egf

#' @rdname getCall.egf
#' @export
#' @importFrom stats getCall
getCall.egf_model_simulate <- function(x, ...) {
  call <- NextMethod("getCall")
  call[[1L]] <- quote(simulate)
  call
}
