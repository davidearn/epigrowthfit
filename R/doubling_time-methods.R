#' Methods for class "doubling_time"
#'
#' @description
#' Methods for "doubling_time" objects returned
#' by [compute_doubling_time()].
#'
#' @param x A "doubling_time" object.
#' @param ... Unused optional arguments.
#'
#' @return
#' The `print` method returns `x` invisibly.
#'
#' @seealso [compute_doubling_time()]
#'
#' @name doubling_time-methods
NULL

#' @rdname doubling_time-methods
#' @export
print.doubling_time <- function(x, ...) {
  cat("Doubling times in days:\n")
  cat("\n")
  print(unclass(x))
  invisible(x)
}
