#' Retrieve calls from objects returned by S3 methods
#'
#' When \code{\link{match.call}} is called from an S3 method for a
#' generic function, a \link{call} to the method is returned, rather
#' than a call to the generic function.
#' Typically, the generic function is exported, but the method is not,
#' so the call cannot be evaluated outside of the namespace environment
#' defining the method.
#' These methods for \code{\link{getCall}} make sure that
#' \code{\link{update.default}} obtains a call to the generic function,
#' rather than the method, when it passes an object to \code{getCall}.
#'
#' @param x
#'   A named \link{list} containing a \link{call} named \code{call}.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' \code{x[["call"]]}, modified so that \code{x[["call"]][[1L]]}
#' is the \link{name} of a generic function, not an S3 method.
#'
#' @examples
#' o1 <- list(call = call("egf.method"))
#' class(o1) <- "egf"
#' getCall(o1)
#' class(o1) <- "egf_no_fit"
#' getCall(o1)
#'
#' o2 <- list(call = call("simulate.method"))
#' class(o2) <- "egf_model_simulate"
#' getCall(o2)
#'
#' @noRd
NULL

#' @export
#' @importFrom stats getCall
getCall.egf <- function(x, ...) {
  call <- NextMethod("getCall")
  call[[1L]] <- quote(egf)
  call
}

#' @export
#' @importFrom stats getCall
getCall.egf_no_fit <- getCall.egf

#' @export
#' @importFrom stats getCall
getCall.egf_model_simulate <- function(x, ...) {
  call <- NextMethod("getCall")
  call[[1L]] <- quote(simulate)
  call
}
