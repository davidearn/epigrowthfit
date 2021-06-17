#' Construct combined model frame
#'
#' Joins in a single \link[=data.frame]{data frame} all mixed effects
#' \link[=model.frame]{model frames} and further variables specified
#' via \code{\link{egf}} argument \code{append}.
#'
#' @param object An \code{"\link{egf}"} object.
#'
#' @details
#' If a variable name occurs in multiple mixed effects model frames,
#' then only one instance is retained. Except in unusual cases
#' (possible only if model formulae have different formula environments),
#' all instances of a variable name are identical, and no information is lost.
#'
#' Since the data frames being combined each correspond rowwise
#' to \code{object$endpoints}, so does the result.
#'
#' @return
#' A \link[=data.frame]{data frame} combining (in the sense
#' of \code{\link{cbind}}) all data frames in the \link{list}
#' \code{object$frame_par} and the data frame \code{object$frame_append}.
#'
#' @keywords internal
make_combined <- function(object) {
  stop_if_not(
    inherits(object, "egf"),
    m = "`object` must inherit from class \"egf\"."
  )
  l <- c(unname(object$frame_par), list(object$frame_append))
  combined <- do.call(cbind, l)
  combined[!duplicated(names(combined))]
}

