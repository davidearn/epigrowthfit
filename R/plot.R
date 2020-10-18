#' Annotate a date axis
#'
#' @description
#' Labels day, month, and year on a horizontal date axis,
#' taking care to ensure that labels are nicely spaced.
#'
#' @param d An increasing Date vector with length 2 or greater.
#'
#' @return
#' A character scalar suggesting a title for the time axis
#' (invisibly).
#'
#' @details
#' The appearance of the axis depends greatly on the range of dates
#' contained in `d`.
#'
#' @keywords internal
#' @export
#' @importFrom graphics axis
daxis <- function(d) {
  if (missing(d)) {
    stop("Missing argument `d`.")
  } else if (!inherits(d, "Date") || length(d) < 2) {
    stop("`d` must be a Date vector with length 2 or greater.")
  } else if (anyNA(d)) {
    stop("`d` must not have missing values.")
  } else if (!all(diff(d) > 0)) {
    stop("`d` must be increasing.")
  }

  r <- as.numeric(diff(range(d)))
  chard <- as.character(d)
  ymd <- matrix(unlist(strsplit(chard, "-")), ncol = 3, byrow = TRUE)
}
