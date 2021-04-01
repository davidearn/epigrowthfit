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
#' #   is.numeric(x),
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
    is.logical(x),
    length(x) == 1L,
    !is.na(x),
    m = sprintf("`%s` must be TRUE or FALSE.", s),
    n = 1L + n
  )
}

stop_if_not_integer <- function(x, n = 1L) {
  s <- deparse(substitute(x))
  stop_if_not(
    is.numeric(x),
    length(x) == 1L,
    x %% 1 == 0,
    m = sprintf("`%s` must be an integer.", s),
    n = 1L + n
  )
}

stop_if_not_positive_integer <- function(x, n = 1L) {
  s <- deparse(substitute(x))
  stop_if_not(
    is.numeric(x),
    length(x) == 1L,
    x >= 1,
    x %% 1 == 0,
    m = sprintf("`%s` must be a positive integer.", s),
    n = 1L + n
  )
}

stop_if_not_number_in_interval <- function(x, a, b, include = c("()", "(]", "[)", "[]"), n = 1L) {
  include <- match.arg(include)
  e1 <- substr(include, 1L, 1L)
  e2 <- substr(include, 2L, 2L)
  f1 <- switch(e1, `(` = `>`, `[` = `>=`)
  f2 <- switch(e2, `)` = `<`, `]` = `<=`)

  s_x <- deparse(substitute(x))
  s_I <- paste0(e1, deparse(substitute(a)), ",",
                deparse(substitute(b)), e2)
  stop_if_not(
    is.numeric(x),
    length(x) == 1L,
    f1(x, a),
    f2(x, b),
    m = sprintf("`%s` must be a number in the interval %s.", s_x, s_I),
    n = 1L + n
  )
}

stop_if_not_character_string <- function(x, n = 1L) {
  s <- deparse(substitute(x))
  stop_if_not(
    is.character(x),
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
#'   For `enum_dupl_string()`, a character vector.
#'   For `enum_dupl_names()`, an \R object with a `names` attribute.
#'
#' @return
#' `enum_dupl_string(x)` returns an enumerated copy of `x`.
#' `enum_dupl_names(x)` returns a copy of `x` with `names(x)`
#' replaced by `enum_dupl_string(names(x))`.
#'
#' @examples
#' # x <- sample(letters[1:3], 10L, replace = TRUE)
#' # enum_dupl_string(x)
#' #
#' # y <- seq_along(x)
#' # names(y) <- x
#' # enum_dupl_names(y)
#'
#' @name enum_dupl_str
#' @keywords internal
NULL

#' @rdname enum_dupl_str
enum_dupl_string <- function(x) {
  f <- factor(x)
  n <- tabulate(f)
  i <- unsplit(lapply(n, seq_len), f)
  sprintf("%s[%d]", x, i)
}

#' @rdname enum_dupl_str
enum_dupl_names <- function(x) {
  `names<-`(x, enum_dupl_string(names(x)))
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
  ok <- !is.na(x)
  x <- x[ok]
  if (length(x) == 0L) {
    return(TRUE)
  }
  if (is.numeric(x)) {
    yes <- abs(max(x) - min(x)) < tol
  } else {
    yes <- length(unique(x)) == 1L
  }
  if (yes && !na.rm && any(!ok)) {
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
    is.numeric(sd),
    length(sd) == d[1L],
    m = "`sd` must be a numeric vector of length `nrow(cor)`."
  )
  cor[] <- sd * cor * rep(sd, each = length(sd))
  cor
}

#' Decompose a vector of dates
#'
#' Extracts year, month, and day from a Date vector.
#'
#' @param x
#'   A Date vector.
#' @param which
#'   A subset of `1:3` indicating which of year, month, and day
#'   should be returned.
#' @param drop
#'   A logical scalar. If `drop = TRUE`, then an integer vector
#'   is returned instead of a matrix when one of `x` and `which`
#'   has length less than 2,
#'
#' @return
#' `X[, which, drop]`, where `X` is an integer matrix with
#' `length(x)` rows and 3 columns listing year, month, and day.
#'
#' @keywords internal
ymd <- function(x, which = 1:3, drop = TRUE) {
  X <- matrix(as.integer(unlist(strsplit(as.character(x), "-"))),
    nrow = length(x),
    ncol = 3L,
    byrow = TRUE,
    dimnames = list(NULL, c("y", "m", "d"))
  )
  X[, which, drop = drop]
}

#' Get ceiling of a date
#'
#' Rounds each Date in a Date vector to the next first-of-the-month
#' or first-of-the-year.
#'
#' @param x A Date vector.
#' @param to A character string.
#'
#' @return
#' `x` with elements replaced by firsts-of-the-month (`"YYYY-MM-01"`)
#' or firsts-of-the-year (`"YYYY-01-01"`), depending on `to`.
#'
#' @keywords internal
dceiling <- function(x, to = c("month", "year")) {
  if (length(x) == 0L) {
    return(.Date(numeric(0L)))
  }
  to <- match.arg(to)
  X <- as.data.frame(ymd(x, drop = FALSE))
  if (to == "month") {
    X$m <- X$m + (X$d > 1L)
    X$y <- X$y + (i <- X$m == 13L)
    X$m[i] <- 1L
    as.Date(paste(X$y, X$m, "1", sep = "-"), format = "%Y-%m-%d")
  } else { # "year"
    X$y <- X$y + (X$m > 1L || X$d > 1L)
    as.Date(paste(X$y, "1", "1", sep = "-"), format = "%Y-%m-%d")
  }
}

#' Recursively merge lists
#'
#' Recursively merge elements of a partially specified list
#' into a completely specified list of defaults.
#'
#' @param x
#'   A list.
#' @param template
#'   A list into which `x` is merged.
#'
#' @details
#' The recursion proceeds as follows:
#'
#' If `x` is `NULL`, then `template` is replaced with `NULL`.
#'
#' Otherwise, if `x` and `template` are both lists or are both
#' atomic vectors and `template` is unnamed, then `x` is recycled
#' to the length of `template`.
#'
#' Otherwise, if `x` and `template` are both lists and `template`
#' is named, then elements of `x` are recursively merged into the
#' so-named elements of `template`.
#'
#' Otherwise, if `x` and `template` are both atomic vectors and
#' `template` is named, then elements of `x` replace the so-named
#' elements of `template`.
#'
#' Otherwise, `x` is discarded and `template` is preserved as is.
#'
#' @return
#' A list.
#'
#' @examples
#' x <- list(
#'   a = NULL,
#'   b = 1,
#'   c = list(
#'     c1 = list(
#'       c11 = 1,
#'       c12 = "a"
#'     ),
#'     c2 = c(Jan = 12L, Dec = 1L),
#'     c3 = FALSE
#'   )
#' )
#' template <- list(
#'   a = list(
#'     a1 = 0,
#'     a2 = "",
#'     a3 = FALSE
#'   ),
#'   b = seq_len(10L),
#'   c = list(
#'     c1 = list(
#'       c11 = 0,
#'       c12 = "",
#'       c13 = FALSE
#'     ),
#'     c2 = `names<-`(seq_len(12L), month.abb),
#'     c3 = list(TRUE)
#'   )
#' )
#'
#' ## Not run
#' \dontrun{
#' rlmerge(x, template)
#' }
#'
#' @keywords internal
rlmerge <- function(x, template) {
  if (is.null(x)) {
    return(NULL)
  }
  if ((is.list(template) && is.list(x)) ||
      (is.atomic(template) && is.atomic(x))) {
    if (is.null(tn <- names(template))) {
      names(x) <- NULL
      return(rep_len(x, length(template)))
    }
    if (!is.null(xn <- names(x))) {
      s <- intersect(tn, xn)
      if (is.list(template)) {
        template[s] <- Map(rlmerge, x = x[s], template = template[s])
      } else {
        template[s] <- x[s]
      }
      return(template)
    }
  }
  template
}

#' Recursively split data frames
#'
#' Splits a supplied data frame on one or more variables in turn.
#'
#' @param x
#'   A data frame.
#' @param by
#'   A subset of `seq_along(x)` or `names(x)` indicating (in order)
#'   variables on which to split `x`.
#' @param drop
#'   A logical scalar passed to [split()].
#'
#' @return
#' A recursive list of data frames with `n = length(by)` levels,
#' unless `n` is zero, in which case `x` is returned as is.
#'
#' @examples
#' f <- function() sample(letters[1:5], 20L, replace = TRUE)
#' d <- as.data.frame(replicate(3L, f()))
#'
#' ## Not run
#' \dontrun{
#' rsplit(d, 1:3)
#' }
#'
#' @keywords internal
rsplit <- function(x, by = integer(0L), drop = FALSE) {
  if (length(by) == 0L) {
    return(x)
  }
  l <- split(x, x[[by[1L]]], drop = drop)
  if (length(by) == 1L) {
    return(l)
  }
  lapply(l, rsplit, by = by[-1L], drop = drop)
}

#' Utilities for nonstandard evaluation
#'
#' Utilities for obtaining index vectors from unevaluated `subset`,
#' `append`, and `order` expressions and character vectors from
#' unevaluated `label` expressions (`xlab`, `ylab`, etc.). Used by
#' several methods for class `"egf"`.
#'
#' @param subset,order,label
#'   Expressions to be evaluated in `data`.
#' @param append
#'   An expression to be evaluated in
#'   ````names<-```(as.list(seq_along(data)), names(data))`.
#' @param data
#'   A data frame.
#' @param enclos
#'   An environment to enclose `data`.
#' @param .subset,.order,.append,`.label`
#'   Atomic vectors to be used (if non-`NULL`) in place of
#'   the result of evaluating the undotted argument.
#'
#' @details
#' `subset` must evaluate to a logical vector of length
#' `n = nrow(data)`.
#' `NULL` is equivalent to `rep_len(TRUE, n)`.
#'
#' `order` must evaluate to a permutation of `seq_len(n)`.
#' `NULL` is equivalent to `seq_len(n)`.
#'
#' `append` must evaluate to an atomic vector indexing `data`.
#' `NULL` is equivalent to `integer(0L)`.
#'
#' `label` must evaluate to an atomic vector of length 1 or `n`.
#' `NULL` is a no-op.
#'
#' Note that `subset` and `append` are processed similarly to
#' [subset()] arguments `subset` and `select`. See [subset()]
#' for additional usage examples.
#'
#' # Warning
#'
#' Nonstandard evaluation of `subset`, `order`, `append`, and `label`
#' is intended only to make interactive use more convenient. To avoid
#' unexpected behaviour, especially when programming, use the dotted
#' versions.
#'
#' @return
#' Let `r` be the result of evaluating the supplied expression
#' or, otherwise, the dotted argument.
#'
#' `subset_to_index()` returns `which(r)`.
#'
#' `order_to_index()` returns `r` as is.
#'
#' `append_to_index()` returns
#' `match(names(data[r]), names(data), 0L)`, without zeros.
#'
#' `label_to_character()` returns
#' `rep_len(as.character(r), nrow(data))`.
#'
#' @examples
#' year <- 2021L
#' data <- data.frame(
#'   month = sample(month.abb, 20L, replace = TRUE),
#'   day = sample(30L, 20L, replace = TRUE)
#' )
#'
#' subset <- substitute(grepl("^J", month) & day < 16L)
#' order <- substitute(order(month, day))
#' append <- substitute(-day)
#' label <- substitute(sprintf("%d-%02d-%02d", year, match(month, month.abb, 0L), day))
#'
#' ## Not run
#' \dontrun{
#' subset_to_index(subset, data, parent.frame())
#' order_to_index(order, data, parent.frame())
#' append_to_index(append, data, parent.frame())
#' label_to_character(label, data, parent.frame())
#' }
#'
#' @name nse
#' @keywords internal
NULL

#' @rdname nse
subset_to_index <- function(subset, data, enclos, .subset = NULL) {
  n <- nrow(data)
  if (is.null(.subset)) {
    if (is.null(subset)) {
      return(seq_len(n))
    }
    r <- eval(subset, data, enclos)
    s <- "`subset` must evaluate to"
  } else {
    r <- .subset
    s <- "`.subset` must be"
  }
  stop_if_not(
    is.logical(r),
    length(r) == n,
    m = paste(s, "a logical vector of length `nrow(data)`.")
  )
  which(r)
}

#' @rdname nse
append_to_index <- function(append, data, enclos, .append = NULL) {
  if (is.null(.append)) {
    if (is.null(append)) {
      return(integer(0L))
    }
    l <- as.list(seq_along(data))
    names(l) <- names(data)
    r <- eval(append, l, enclos)
  } else {
    r <- .append
  }
  m <- match(names(data[r]), names(data), 0L)
  m[m > 0L]
}

#' @rdname nse
order_to_index <- function(order, data, enclos, .order = NULL) {
  n <- nrow(data)
  if (is.null(.order)) {
    if (is.null(order)) {
      return(seq_len(n))
    }
    r <- eval(order, data, enclos)
    s <- "`order` must evaluate to "
  } else {
    r <- .order
    s <- "`.order` must be "
  }
  stop_if_not(
    is.numeric(r),
    length(r) == n,
    sort(r) == seq_len(n),
    m = paste0(s, "a permutation\nof `seq_len(nrow(data))`.")
  )
  r
}

#' @rdname nse
label_to_character <- function(label, data, enclos, .label = NULL) {
  n <- nrow(data)
  if (is.null(.label)) {
    if (is.null(label)) {
      return(NULL)
    }
    r <- eval(label, data, enclos)
    s <- "`label` must evaluate to "
  } else {
    r <- .label
    s <- "`.label` must be "
  }
  stop_if_not(
    is.atomic(label),
    any(length(label) == c(1L, n)),
    m = paste0(s, "an atomic vector\nof length 1 or `nrow(data)`.")
  )
  rep_len(as.character(r), n)
}
