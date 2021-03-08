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

stop_if_not_true_false <- function(x, n = 1L) {
  s <- deparse(substitute(x))
  stop_if_not(
    is.vector(x, "logical"),
    length(x) == 1L,
    !is.na(x),
    m = sprintf("`%s` must be TRUE or FALSE.", s),
    n = 1L + n
  )
}

stop_if_not_integer <- function(x, n = 1L) {
  s <- deparse(substitute(x))
  stop_if_not(
    is.vector(x, "numeric"),
    length(x) == 1L,
    x %% 1 == 0,
    m = sprintf("`%s` must be an integer.", s),
    n = 1L + n
  )
}

stop_if_not_positive_integer <- function(x, n = 1L) {
  s <- deparse(substitute(x))
  stop_if_not(
    is.vector(x, "numeric"),
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
    is.vector(x, "numeric"),
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
    is.vector(x, "character"),
    length(x) == 1L,
    !is.na(x),
    m = sprintf("`%s` must be a character string.", s),
    n = 1L + n
  )
}

#' Enumerate duplicated strings
#'
#' Replaces the `i`th instance of string `s` in character vector `x`
#' with `sprintf("%s[%d]", s, i)`.
#'
#' @param x A character vector.
#'
#' @return
#' `x` with elements enumerated.
#'
#' @examples
#' x <- sample(letters[1:3], 10L, replace = TRUE)
#' enum_dupl_str(x)
#'
#' @noRd
enum_dupl_str <- function(x) {
  f <- factor(x)
  n <- tabulate(f)
  i <- unsplit(lapply(n, seq_len), f)
  sprintf("%s[%d]", x, i)
}



#' Test whether a vector is "constant"
#'
#'
#'
is_constant <- function(x, na.rm = FALSE, tol = sqrt(.Machine$double.eps)) {
  if (is.data.frame(x)) {
    return(all(vapply(x, is_constant, FALSE, na.rm, tol)))
  }
  if (length(x) == 0L || (na.rm && all(is.na(x)))) {
    return(TRUE)
  }
  if (is.numeric(x)) {
    yes <- abs(max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) < tol
  } else {
    yes <- length(unique(x[!is.na(x)])) == 1L
  }
  if (yes && anyNA(x) && !na.rm) {
    return(NA)
  }
  yes
}
