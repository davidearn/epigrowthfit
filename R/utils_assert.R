#' Utilities for object validation
#'
#' @param ...
#'   \link[=logical]{Logical} vectors.
#' @param m
#'   A \link{character} string used as an error or warning message.
#' @param n
#'   A non-negative integer specifying a call to be associated with
#'   the error or warning message. \code{n} indicates a number of
#'   generations backwards relative to the function call. The default
#'   (\code{n = 1}) corresponds to the call from the parent frame.
#'
#' @details
#' \code{stop_if_not} is a replacement for \code{\link{stopifnot}}
#' allowing the user to specify an error message and function call.
#'
#' \code{warn_if_not} behaves identically to \code{stop_if_not} but
#' issues a warning instead of an error.
#'
#' @return
#' \code{\link{NULL}} (invisibly).
#'
#' @examples
#' ## x <- 1.1
#' ## stop_if_not(
#' ##   is.numeric(x),
#' ##   length(x) == 1L,
#' ##   x > 0,
#' ##   x < 1,
#' ##   m = "`x` must be a number in the interval (0,1)."
#' ## )
#'
#' @name stop_if_not
#' @keywords internal
NULL

#' @rdname stop_if_not
stop_if_not <- function(..., m = "", n = 1L) {
  l <- ...length()
  for (i in seq_len(l)) {
    e <- ...elt(i)
    if (!(is.logical(e) && !anyNA(e) && all(e))) {
      p <- sys.parent(n)
      stop(simpleError(m, call = if (p > 0L) sys.call(p)))
    }
  }
  invisible(NULL)
}

#' @rdname stop_if_not
warn_if_not <- function(..., m = "", n = 1L) {
  l <- ...length()
  for (i in seq_len(l)) {
    e <- ...elt(i)
    if (!(is.logical(e) && !anyNA(e) && all(e))) {
      p <- sys.parent(n)
      warning(simpleWarning(m, call = if (p > 0L) sys.call(p)))
    }
  }
  invisible(NULL)
}

stop_if_not_true_false <- function(x, allow_numeric = FALSE, n = 1L) {
  s <- deparse(substitute(x))
  a <- if (allow_numeric) " or a number" else ""
  stop_if_not(
    is.logical(x) || (allow_numeric && is.numeric(x)),
    length(x) == 1L,
    !is.na(x),
    m = sprintf("`%s` must be TRUE or FALSE%s.", s, a),
    n = 1L + n
  )
}

stop_if_not_integer <- function(x, kind = c("any", "positive", "nonnegative", "negative", "nonpositive"), n = 1L) {
  s <- deparse(substitute(x))
  kind <- match.arg(kind)
  a <- switch(kind, any = "an", paste("a", kind))
  f <- switch(kind,
    any         = function(x, y) TRUE,
    positive    = `>`,
    nonnegative = `>=`,
    negative    = `<`,
    nonpositive = `<=`
  )
  stop_if_not(
    is.numeric(x),
    length(x) == 1L,
    x %% 1 == 0,
    f(x, 0),
    m = sprintf("`%s` must be %s integer.", s, a),
    n = 1L + n
  )
}

stop_if_not_number <- function(x, kind = c("any", "positive", "nonnegative", "negative", "nonpositive"), n = 1L) {
  s <- deparse(substitute(x))
  kind <- match.arg(kind)
  a <- switch(kind, any = "a", paste("a", kind))
  f <- switch(kind,
    any         = function(x, y) TRUE,
    positive    = `>`,
    nonnegative = `>=`,
    negative    = `<`,
    nonpositive = `<=`
  )
  stop_if_not(
    is.numeric(x),
    length(x) == 1L,
    f(x, 0),
    m = sprintf("`%s` must be %s number.", s, a),
    n = 1L + n
  )
}

stop_if_not_number_in_interval <- function(x, a = -Inf, b = Inf, include = c("()", "(]", "[)", "[]"), n = 1L) {
  s <- deparse(substitute(x))
  include <- match.arg(include)
  d1 <- substr(include, 1L, 1L)
  d2 <- substr(include, 2L, 2L)
  interval <- paste0(d1, deparse(substitute(a)), ",", deparse(substitute(b)), d2)
  f1 <- switch(d1, `(` = `>`, `[` = `>=`)
  f2 <- switch(d2, `)` = `<`, `]` = `<=`)
  stop_if_not(
    is.numeric(x),
    length(x) == 1L,
    f1(x, a),
    f2(x, b),
    m = sprintf("`%s` must be a number in the interval %s.", s, interval),
    n = 1L + n
  )
}

stop_if_not_string <- function(x, n = 1L) {
  s <- deparse(substitute(x))
  stop_if_not(
    is.character(x),
    length(x) == 1L,
    !is.na(x),
    m = sprintf("`%s` must be a character string.", s),
    n = 1L + n
  )
}

stop_if_not_Date <- function(x, n = 1L) {
  s <- deparse(substitute(x))
  stop_if_not(
    inherits(x, "Date"),
    length(x) == 1L,
    !is.na(x),
    m = sprintf("`%s` must be a Date vector of length one.", s),
    n = 1L + n
  )
}
