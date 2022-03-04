#' Utilities for object validation
#'
#' \code{stop_if_not} is a drop-in replacement for \code{\link{stopifnot}},
#' signalling an error with a user-specified message if at least one assertion
#' is not true. \code{warn_if_not} behaves identically but issues a warning
#' instead of an error.
#'
#' @param ...
#'   Logical vectors.
#' @param m
#'   A \link{character} string used as an error or warning message.
#' @param n
#'   A non-negative integer specifying a call to be associated with the error
#'   or warning message. \code{n} indicates a number of generations backwards
#'   relative to the evaluation frame.
#'   The default (\code{n = 1}) corresponds to the call from the parent frame.
#'
#' @return
#' \code{NULL} (invisibly).
#'
#' @examples
#' x <- 1.1
#' stop_if_not(is.double(x), length(x) == 1L, x > 0, x < 1,
#'             m = "'x' must be a number in the interval (0,1).")
#'
#' @noRd
NULL

stopifnot1 <- function(..., m = "", n = 1L) {
    for (i in seq_len(...length())) {
        e <- ...elt(i)
        if (!(is.logical(e) && !anyNA(e) && all(e))) {
            p <- sys.parent(n)
            stop(simpleError(m, call = if (p > 0L) sys.call(p)))
        }
    }
    invisible(NULL)
}

warnifnot1 <- function(..., m = "", n = 1L) {
    for (i in seq_len(...length())) {
        e <- ...elt(i)
        if (!(is.logical(e) && !anyNA(e) && all(e))) {
            p <- sys.parent(n)
            warning(simpleWarning(m, call = if (p > 0L) sys.call(p)))
        }
    }
    invisible(NULL)
}

stop1 <- function(...) {
    m <- wrap(...)
    p <- sys.parent(1L)
    stop(simpleError(m, call = if (p > 0L) sys.call(p)))
}

warning1 <- function(...) {
    m <- wrap(...)
    p <- sys.parent(1L)
    warning(simpleWarning(m, call = if (p > 0L) sys.call(p)))
}

is_true_or_false <- function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x)
}

is_flag <- function(x) {
    (is.logical(x) || is.numeric(x)) && length(x) == 1L && is.finite(x)
}

is_number <- function(x,
                      kind = c("any",
                               "positive", "nonnegative",
                               "negative", "nonpositive"),
                      integer = FALSE) {
    kind <- match.arg(kind)
    relop <- switch(kind,
                    any = function(x, y) TRUE,
                    positive = `>`,
                    nonnegative = `>=`,
                    negative = `<`,
                    nonpositive = `<=`)
    is.numeric(x) && length(x) == 1L && is.finite(x) &&
        relop(x, 0) && (!integer || x %% 1 == 0)
}

is_number_in_interval <- function(x, a = -Inf, b = Inf,
                                  include = c("()", "(]", "[)", "[]")) {
    include <- match.arg(include)
    relop1 <- switch(substr(include, 1L, 1L), `(` = `>`, `[` = `>=`)
    relop2 <- switch(substr(include, 2L, 2L), `)` = `<`, `]` = `<=`)
    is.numeric(x) && length(x) == 1L && !is.na(x) &&
        relop1(x, a) && relop2(x, b)
}

is_string <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x)
}
