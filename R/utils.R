#' Validate an atomic vector
#'
#' @description
#' Performs checks on objects (usually vectors) and optionally
#' throws an error if desired conditions are not met. Simplifies
#' the task of validating atomic vector arguments of functions.
#'
#' @param x
#'   An object.
#' @param ...
#'   Zero or more objects coercible to character, which are
#'   pasted together without separator to generate an error
#'   or warning message (depending on `action`).
#' @param what
#'   A character vector listing object classes. Check fails
#'   if `x` does not inherit from one of these classes.
#' @param not
#'   A character vector listing object classes. Check fails
#'   if `x` inherits from one of these classes.
#' @param len
#'   A numeric vector. If `length(len) = 1`, then check fails
#'   if `length(x) != len`. If `length(len) > 1`, then check
#'   fails if `!is.na(len[1])` and `length(x) < len[1]`,
#'   or if `!is.na(len[2])` and `length(x) > len[2]`.
#' @param opt
#'   A vector giving a discrete set of acceptable values
#'   for `x[[i]]`. Check fails if `!all(x %in% opt)`.
#' @param val
#'   A numeric vector specifying an interval of acceptable
#'   values for `x[i]`, assuming that `x` is numeric.
#'   If `length(val) = 1`, then check fails if
#'   `!all(x == val, na.rm = TRUE)`.
#'   If `length(val) > 1`, then `val[1]` and `val[2]` are
#'   interpreted as interval endpoints, and whether the
#'   endpoints are included in the interval is determined
#'   by `rel`. For example, if `rel = c(">=", "<=")`,
#'   then check fails if
#'   `!is.na(val[1])` and `!all(x >= val[1], na.rm = TRUE)`,
#'   or if `!is.na(val[2])` and `!all(x <= val[2], na.rm = TRUE)`.
#' @param rel
#'   A character vector of length 2 such that
#'   `rel[1] %in% c(">", ">=")` and `rel[2] %in% c("<", "<=")`.
#' @param yes
#'   A function or list of functions of one argument.
#'   Check fails if any does not return `TRUE` when applied to `x`.
#' @param no
#'   A function or list of functions of one argument.
#'   Check fails if any does not return `FALSE` when applied to `x`.
#' @param action
#'   One of `"stop"`, `"warning"`, and `"none"`,
#'   indicating how check failure should be handled.
#'
#' @return
#' `TRUE` if check passes and `FALSE` otherwise.
#'
#' @details
#' Checks are performed in the order of `names(formals(check))`.
#' If `what` fails then `not` will not be checked. If `not` fails,
#' then `len` won't be checked (and so on).
#'
#' Ironically, no checks on the input are performed, so use with
#' care (e.g., if setting `val`, then make sure that `x` is numeric
#' or that `what` contains `"numeric"`).
#'
#' Integer vectors are accepted as numeric vectors,
#' hence `check(0L, what = "numeric")` returns `TRUE`.
#'
#' @keywords internal
check <- function(x, ..., what = NULL, not = NULL, len = NULL,
                  opt = NULL, val = NULL, rel = c(">=", "<="),
                  yes = NULL, no = NULL,
                  action = c("stop", "warning", "none")) {
  action <- match.arg(action)

  ## Tolerate integer when checking for numeric
  if ("numeric" %in% what) {
    what <- c(what, "integer")
  }
  ## Tolerate function instead of list of functions
  if (is.function(yes)) {
    yes <- list(yes)
  }
  if (is.function(no)) {
    no <- list(no)
  }

  passes <-
    ## Check class
    (is.null(what) || inherits(x, what)) &&
    (is.null(not) || !inherits(x, not)) &&
    ## Check length: equals `len` or in `range(len)`
    (is.null(len) ||
       (length(len) == 1L && length(x) == len) ||
       (length(len) > 1L &&
          (is.na(len[1L]) || length(x) >= len[1L]) &&
          (is.na(len[2L]) || length(x) <= len[2L]))) &&
    ## Check value: is a subset of `opt`
    (is.null(opt) || all(x %in% opt)) &&
    ## Check value: equals `val` or in `range(val)`
    (is.null(val) ||
       (length(val) == 1L && all(x == val, na.rm = TRUE)) ||
       (length(val) > 1L &&
          (is.na(val[1L]) || all(match.fun(rel[1L])(x, val[1L]), na.rm = TRUE)) &&
          (is.na(val[2L]) || all(match.fun(rel[2L])(x, val[2L]), na.rm = TRUE)))) &&
    ## Check value: passes functions `yes`
    (is.null(yes) || all(sapply(yes, function(f) isTRUE(f(x))))) &&
    ## Check value: fails functions `no`
    (is.null(no) || all(sapply(no, function(f) isFALSE(f(x)))))
  if (!isTRUE(passes) && action != "none") {
    do.call(action, list(..., call. = FALSE))
  }
  passes
}
