#' Utilities for object validation
#'
#' @param ...
#'   Logical vectors.
#' @param m
#'   A character string used as an error or warning message.
#' @param n
#'   A non-negative integer specifying a call to be associated with
#'   the error or warning. `n` indicates a number of generations
#'   backwards relative to the function call. The default `n = 1L`
#'   corresponds to the call from the parent frame.
#' @param x
#'   An object to be checked.
#'
#' @details
#' `stop_if_not()` is a replacement for [base::stopifnot()] allowing
#' (requiring) the user to specify an error message and function call.
#'
#' `warn_if_not()` behaves identically to `stop_if_not()` but issues
#' a warning instead of an error.
#'
#' `stop_if_not_tf()` stops if its argument is not `TRUE` or `FALSE`.
#' `stop_if_not_positive_integer()` stops if its argument is not a
#' positive integer (of type `"integer"` or `"double"`).
#' `stop_if_not_in_0_1()` stops if its arguments is not a number in
#' the interval (0,1).
#' `stop_if_not_character_string()` stops if its argument is not a
#' character vector of length 1 (other than `NA_character_`).
#'
#' @return
#' `NULL` (invisibly).
#'
#' @noRd
NULL

stop_if_not <- function(..., m = "", n = 1L) {
  n <- ...length()
  for (i in seq_len(n)) {
    e <- ...elt(i)
    if (!(is.logical(e) && !anyNA(e) && all(e))) {
      p <- sys.parent(n)
      stop(simpleError(m, call = if (p > 0L) sys.call(p)))
    }
  }
  invisible(NULL)
}

warn_if_not <- function(..., m, n = 1L) {
  n <- ...length()
  for (i in seq_len(n)) {
    e <- ...elt(i)
    if (!(is.logical(e) && !anyNA(e) && all(e))) {
      p <- sys.parent(n)
      warning(simpleWarning(m, call = if (p > 0L) sys.call(p)))
    }
  }
  invisible(NULL)
}

stop_if_not_tf <- function(x, n = 1L) {
  s <- deparse(substitute(x))
  stop_if_not(
    inherits(x, "logical"),
    length(x) == 1L,
    !is.na(x),
    m = sprintf("`%s` must be TRUE or FALSE.", s),
    n = 1L + n
  )
}

stop_if_not_positive_integer <- function(x, n = 1L) {
  s <- deparse(substitute(x))
  stop_if_not(
    inherits(x, c("integer", "numeric")),
    length(x) == 1L,
    x >= 1,
    x %% 1 == 0,
    m = sprintf("`%s` must be a positive integer.", s),
    n = 1L + n
  )
}

stop_if_not_in_0_1 <- function(x, n = 1L) {
  s <- deparse(substitute(x))
  stop_if_not(
    inherits(x, "numeric"),
    length(x) == 1L,
    x > 0,
    x < 1,
    m = sprintf("`%s` must be a number in the interval (0,1).", s),
    n = 1L + n
  )
}

stop_if_not_character_string <- function(x, n = 1L) {
  s <- deparse(substitute(x))
  stop_if_not(
    inherits(x, "character"),
    length(x) == 1L,
    !is.na(x),
    m = sprintf("`%s` must be a character string.", s),
    n = 1L + n
  )
}

is_constant_within_tol <- function(x, tol = sqrt(.Machine$double.eps), na.rm = FALSE) {
  if (length(x) == 0L || (na.rm && all(is.na(x)))) {
    return(NA)
  }
  abs(max(x, na.rm = na.rm) - min(x, na.rm = na.rm)) < tol
}

is_integer_within_tol <- function(x, tol = sqrt(.Machine$double.eps)) {
  abs(x - round(x)) < tol
}
