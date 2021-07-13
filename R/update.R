#' Update a fitted model
#'
#' When \code{\link{match.call}} is called from a method for a generic function,
#' a \link{call} to the method is returned, rather than a call to the
#' generic function. Typically, the method is not exported, so the call
#' cannot be evaluated by the default method for \code{\link{update}}.
#' These methods for \code{\link{update}} circumvent this issue for objects
#' returned by \code{\link{egf}}.
#'
#' @param object
#'   An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.
#' @param ...
#'   Arguments for the updated \link{call} (see \code{\link{update}}).
#' @param evaluate
#'   A \link{logical} flag, indicating whether the updated \link{call}
#'   should be evaluated.
#'
#' @return
#' The updated \link{call} or the result of its evaluation,
#' depending on \code{evaluate}.
#'
#' @name update.egf
NULL

#' @rdname update.egf
#' @export
update.egf <- function(object, ..., evaluate = TRUE) {
  object$call[[1L]] <- quote(egf)
  NextMethod("update")
}

#' @rdname update.egf
#' @export
update.egf_no_fit <- update.egf
