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
#' `stop_if_not()` is a replacement for [stopifnot()] allowing
#' (requiring) the user to specify an error message and function call.
#'
#' `warn_if_not()` behaves identically to `stop_if_not()` but issues
#' a warning instead of an error.
#'
#' @return
#' `NULL` (invisibly).
#'
#' @examples
#' # x <- 1.1
#' # stop_if_not(
#' #   is.vector(x, "numeric"),
#' #   length(x) == 1L,
#' #   x > 0,
#' #   x < 1,
#' #   m = "`x` must be a number in the interval (0,1)."
#' # )
#'
#' @name stop_if_not
#' @keywords internal
NULL

#' @rdname stop_if_not
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

#' @rdname stop_if_not
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

stop_if_not_number_in_interval <- function(x, a, b, include = c("()", "(]", "[)", "[]"), n = 1L) {
  s <- deparse(substitute(x))
  include <- match.arg(include)
  e1 <- substr(include, 1L, 1L)
  e2 <- substr(include, 2L, 2L)
  f1 <- switch(e1, `(` = `>`, `[` = `>=`)
  f2 <- switch(e2, `)` = `<`, `]` = `<=`)
  interval <- paste0(e1, deparse(substitute(a)), ",", deparse(substitute(b)), e2)
  stop_if_not(
    is.vector(x, "numeric"),
    length(x) == 1L,
    f1(x, a),
    f2(x, b),
    m = sprintf("`%s` must be a number in the interval %s.", s, interval),
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
#' Replaces the `i`th instance of strings `s` in a character vector
#' with `sprintf("%s[%d]", s, i)`.
#'
#' @param x
#'   For `enum_dupl_str()`, a character vector.
#'   For `enum_dupl_names()`, an \R object with a `names` attribute.
#'
#' @return
#' `enum_dupl_str(x)` returns an enumerated copy of `x`.
#' `enum_dupl_names(x)` returns a copy of `x` with `names(x)`
#' replaced by `enum_dupl_str(names(x))`.
#'
#' @examples
#' # x <- sample(letters[1:3], 10L, replace = TRUE)
#' # enum_dupl_str(x)
#' #
#' # y <- seq_along(x)
#' # names(y) <- x
#' # enum_dupl_names(y)
#'
#' @name enum_dupl_str
#' @keywords internal
NULL

#' @rdname enum_dupl_str
enum_dupl_str <- function(x) {
  f <- factor(x)
  n <- tabulate(f)
  i <- unsplit(lapply(n, seq_len), f)
  sprintf("%s[%d]", x, i)
}

#' @rdname enum_dupl_str
enum_dupl_names <- function(x) {
  `names<-`(x, enum_dupl_str(names(x)))
}

#' Test whether an atomic vector is "constant"
#'
#' Tests whether the elements of an atomic vector are equal
#' (or nearly equal in the case of numeric vectors).
#'
#' @param x
#'   An atomic vector or data frame. For data frames `x`,
#'   `is_constant(x)` tests whether all listed vectors are
#'   (individually) constant.
#' @param na.rm
#'   A logical scalar. If `TRUE`, then `NA` in vectors `x` are ignored.
#' @param tol
#'   A positive number. Numeric vectors `x` are considered constant
#'   if and only if the distance from `min(x)` to `max(x)` is less
#'   than `tol`.
#'
#' @details
#' `TRUE` is returned for zero-length `x`.
#'
#' @return
#' `TRUE`, `FALSE`, or `NA`.
#'
#' @examples
#' # x <- c(0, 1e-03, NA)
#' # is_constant(x, na.rm = TRUE, tol = 1e-02)
#' # is_constant(x, na.rm = TRUE, tol = 1e-04)
#' # is_constant(x, na.rm = FALSE, tol = 1e-02)
#' # is_constant(x, na.rm = FALSE, tol = 1e-04)
#'
#' @keywords internal
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
  if (yes && !na.rm && anyNA(x)) {
    return(NA)
  }
  yes
}

#' Correlation to covariance matrix conversion
#'
#' Perform the inverse of [stats::cov2cor()].
#'
#' @param cor
#'   A symmetric numeric matrix. Diagonal elements should be 1 and
#'   off-diagonal elements should be numbers in the interval `[-1,1]`.
#' @param sd
#'   A numeric vector of length `nrow(cor)`.
#'
#' @return
#' The result of `sweep(sweep(cor, 1L, sd, ```*```), 2L, sd, ```*```)`,
#' obtained slightly more efficiently, following [stats::cov2cor()].
#'
#' @examples
#' # X <- replicate(6L, rnorm(10L))
#' # V <- cov(X, X)
#' # all.equal(V, cor2cov(cov2cor(V), sqrt(diag(V))))
#'
#' @keywords internal
cor2cov <- function(cor, sd) {
  d <- dim(cor)
  stop_if_not(
    length(d) == 2L,
    is.numeric(cor),
    d[1L] == d[2L],
    m = "`cor` must be a square numeric matrix."
  )
  stop_if_not(
    is.vector(sd, "numeric"),
    length(sd) == d[1L],
    m = "`sd` must be a numeric vector of length `nrow(cor)`."
  )
  cor[] <- sd * cor * rep(sd, each = length(sd))
  cor
}
