#' Update an object returned by an S3 method
#'
#' When \code{\link{match.call}} is called from an S3 method for a generic
#' function, a \link{call} to the method is returned, rather than a call
#' to the generic function. Typically, the method is not exported,
#' so the call cannot be evaluated by \code{\link{update.default}}.
#' These methods for \code{\link{getCall}} ensure that
#' \code{\link{update.default}} obtains a call to the generic function
#' that it \emph{can} evaluate.
#'
#' @param x
#'   A \link{list} with an element \code{call}, which should be a \link{call}.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' \code{x[["call"]]}, modified so that \code{x[["call"]][[1L]]}
#' is the \link{name} of a generic function, not an S3 method.
#'
#' @examples
#' object <- list(call = call("generic.method"))
#' class(object) <- "egf"
#' getCall(object)
#' class(object) <- "egf_no_fit"
#' getCall(object)
#' class(object) <- "egf_model_simulate"
#' getCall(object)
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
