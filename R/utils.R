#' Enumerate duplicated strings
#'
#' Replaces the \code{i}th instance of string \code{s} in a
#' \link{character} vector with \code{\link{sprintf}("\%s[\%d]", s, i)}.
#'
#' @param x
#'   For \code{enum_dupl_string}, a \link{character} vector.
#'   For \code{enum_dupl_names}, an \R object with a \code{\link{names}}
#'   \link[=attributes]{attribute}.
#'
#' @return
#' \code{enum_dupl_string(x)} returns an enumerated copy of \code{x}.
#' \code{enum_dupl_names(x)} returns a copy of \code{x} with
#' \code{\link{names}(x)} replaced by \code{enum_dupl_string(\link{names}(x))}.
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
#' Tests whether the elements of an \link{atomic} vector are equal
#' (or nearly equal in the case of \link{double} vectors).
#'
#' @param x
#'   An \link{atomic} vector or \link[=data.frame]{data frame}.
#' @param na.rm
#'   A \link{logical} scalar. If \code{TRUE}, then \code{\link{NA}}
#'   are ignored.
#' @param tol
#'   A positive number. \link[=double]{Double} vectors \code{x}
#'   are considered constant if and only if the distance from
#'   \code{\link{min}(x)} to \code{\link{max}(x)} is less than
#'   \code{tol}.
#'
#' @details
#' \code{TRUE} is returned if \code{\link{length}(x) = 0}.
#'
#' For data frames \code{x}, \code{is_constant(x)} tests
#' whether all listed vectors are constant.
#'
#' @return
#' \code{TRUE}, \code{FALSE}, or \code{\link{NA}}.
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
    return(all(vapply(x, is_constant, FALSE, na.rm = na.rm, tol = tol)))
  }
  ok <- !is.na(x)
  x <- x[ok]
  if (length(x) == 0L) {
    return(TRUE)
  }
  if (is.double(x)) {
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
#' Perform the inverse of \code{\link{cov2cor}}.
#'
#' @param cor
#'   A square \link{numeric} \link{matrix}.
#' @param sd
#'   A \link{numeric} vector of length \code{\link{nrow}(cor)}.
#'
#' @return
#' The result of
#' \code{\link{sweep}(\link{sweep}(cor, 1L, sd, `*`), 2L, sd, `*`)},
#' obtained slightly more efficiently, following \code{\link{cov2cor}}.
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

#' Concatenate, wrap, collapse
#'
#' Concatenates \R objects, formats the resulting string in
#' a paragraph using \code{\link{strwrap}}, and collapses
#' paragraph lines with newlines. Under default settings,
#' the result is a string that prints in the console as a
#' nicely wrapped paragraph.
#'
#' @param ...
#'   Zero or more \R objects, to be coerced to \link{character}
#'   and concatenated with no separator.
#' @param width
#'   A positive integer passed to \code{\link{strwrap}},
#'   indicating a target column width.
#'
#' @return
#' A \link{character} vector of length 1.
#'
#' @keywords internal
wrap <- function(..., width = 0.9 * getOption("width")) {
  x <- paste0(..., collapse = " ")
  paragraph_lines <- strwrap(x, width = width)
  paste(paragraph_lines, collapse = "\n")
}

rle_patch <- function(x) {
  n <- length(x)
  if (n == 0L) {
    return(list(lengths = integer(0L), values = x))
  }
  l <- x[-n] != x[-1L]
  if (any(argna <- is.na(x))) {
    l[is.na(l)] <- FALSE
    l <- l | (argna[-n] & !argna[-1L]) | (!argna[-n] & argna[-1L])
  }
  i <- c(which(l), n)
  list(lengths = diff(c(0L, i)), values = x[i])
}

locf <- function(x, x0 = NULL) {
  if (!anyNA(x)) {
    return(x)
  }
  rle_x <- rle_patch(x)
  y <- rle_x$values
  if (is.na(y[1L]) && !is.null(x0)) {
    y[1L] <- x0
  }
  if (anyNA(y[-1L])) {
    argna_y <- which(c(FALSE, is.na(y[-1L])))
    y[argna_y] <- y[argna_y - 1L]
  }
  rle_x$values <- y
  inverse.rle(rle_x)
}

mftapply <- function(x, index, f) {
  stop_if_not(
    is.list(f),
    length(f) > 0L,
    vapply(f, is.function, FALSE),
    m = "`f` must be a list of functions of positive length."
  )
  if (is.data.frame(x)) {
    do_call <- function(f0, x) { x[] <- lapply(x, f0); x }
  } else {
    do_call <- function(f0, x) { f0(x) }
  }
  split(x, index, drop = FALSE) <- Map(do_call, x = split(x, index, drop = FALSE), f0 = f)
  x
}

#' @importFrom stats qchisq
do_wald <- function(estimate, se, level) {
  q <- qchisq(level, df = 1)
  n <- length(estimate)
  lu <- estimate + rep.int(sqrt(q) * c(-1, 1), c(n, n)) * se
  dim(lu) <- c(n, 2L)
  colnames(lu) <- c("lower", "upper")
  lu
}

