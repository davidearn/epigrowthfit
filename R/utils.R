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
#' ## x <- sample(letters[1:3], 10L, replace = TRUE)
#' ## enum_dupl_string(x)
#' ##
#' ## y <- seq_along(x)
#' ## names(y) <- x
#' ## enum_dupl_names(y)
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
#' ## x <- c(0, 1e-03, NA)
#' ## is_constant(x, na.rm = TRUE, tol = 1e-02)
#' ## is_constant(x, na.rm = TRUE, tol = 1e-04)
#' ## is_constant(x, na.rm = FALSE, tol = 1e-02)
#' ## is_constant(x, na.rm = FALSE, tol = 1e-04)
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
#'   A symmetric \link{numeric} \link{matrix}.
#' @param sd
#'   A \link{numeric} vector of length \code{\link{nrow}(cor)}.
#'
#' @return
#' The result of
#' \code{\link{sweep}(\link{sweep}(cor, 1L, sd, `*`), 2L, sd, `*`)},
#' obtained slightly more efficiently, following \code{\link{cov2cor}}.
#'
#' @examples
#' ## X <- replicate(6L, rnorm(10L))
#' ## V <- cov(X, X)
#' ## all.equal(V, cor2cov(cov2cor(V), sqrt(diag(V))))
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
#' Concatenates \R objects, formats the resulting string
#' in paragraphs using \code{\link{strwrap}}, and collapses
#' paragraph lines with newlines. Under default settings,
#' the result is a string that prints in the console as
#' one or more nicely wrapped paragraphs.
#'
#' @param ...
#'   Zero or more \R objects (typically atomic vectors of length 1),
#'   to be coerced to \link{character} and concatenated with no separator.
#' @param width
#'   A positive integer passed to \code{\link{strwrap}},
#'   indicating a target column width.
#'
#' @return
#' A \link{character} vector of length 1.
#'
#' @keywords internal
wrap <- function(..., width = 0.9 * getOption("width")) {
  dots <- unname(list(...))
  x <- vapply(dots, paste0, "", collapse = "")
  x1 <- paste0(x, collapse = "")
  y <- unlist(strsplit(x1, "\n[ \t\n]*\n"), FALSE, FALSE)
  y1 <- paste0(y, collapse = "\n\n")
  z <- strwrap(y1, width = width)
  paste0(z, collapse = "\n")
}

#' Literal run length encoding
#'
#' A replacement for \code{\link{rle}} that regards \code{\link{NA}}
#' as equal to previous \code{NA}. For \link{double} vectors,
#' it regards \code{\link{NaN}} as equal to previous \code{NaN}
#' but unequal to previous \code{NA}.
#'
#' @param x An \link{atomic} vector.
#'
#' @details
#' The result is compatible with \code{\link{inverse.rle}}.
#' More precisely,
#' \code{function(x) inverse.rle(rle_literal(x))}
#' is an identity function for atomic \code{x}, and
#' \code{function(x) rle_literal(inverse.rle(x))}
#' is an identity function for lists \code{x} generated
#' by \code{rle_literal}.
#'
#' @return
#' A \link{list} similar to \code{\link{rle}(x)},
#' though not inheriting from \link{class} \code{"rle"}.
#'
#' @examples
#' ## x <- rep.int(c(0, NA, NaN, 1), 1:4)
#' ## rle_x <- rle_literal(x)
#' ## identical(x, inverse.rle(rle_x))
#'
#' @keywords internal
rle_literal <- function(x) {
  n <- length(x)
  if (n == 0L) {
    return(list(lengths = integer(0L), values = x))
  }
  l <- x[-n] != x[-1L]
  if (any(argna <- is.na(x))) {
    l[is.na(l)] <- FALSE
    if (is.double(x)) {
      argnan <- is.nan(x)
      argna <- argna & !argnan
      l <- l | (argnan[-n] & !argnan[-1L]) | (!argnan[-n] & argnan[-1L])
    }
    l <- l | (argna[-n] & !argna[-1L]) | (!argna[-n] & argna[-1L])
  }
  i <- c(which(l), n)
  list(lengths = diff(c(0L, i)), values = x[i])
}

#' Last observation carried forward
#'
#' Replaces missing values (\code{\link{NA}} and, in \link{double}
#' vectors, \code{\link{NaN}}) with most recent earlier observations.
#'
#' @param x
#'   An \link{atomic} vector.
#' @param x0
#'   An \link{atomic} scalar to replace leading missing values,
#'   for which there are no earlier observations. The default
#'   (\code{\link{NULL}}) is to preserve these missing values.
#'
#' @return
#' \code{x} with missing values replaced.
#'
#' @examples
#' ## x <- c(NA, NA, 1, NA, 2, 2, 3, NA)
#' ## locf(x)
#' ## locf(x, x0 = 0)
#'
#' @keywords internal
locf <- function(x, x0 = NULL) {
  if (!anyNA(x)) {
    return(x)
  }
  rle_x <- rle_literal(x)
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

#' Apply length-preserving functions
#'
#' Replaces subsets of ragged vectors and data frames with the results
#' of length-preserving functions.
#'
#' @param x
#'   A \link{vector} or \link[=data.frame]{data frame}.
#' @param index
#'   A \link{factor} defining a partition of the elements or rows of \code{x}.
#'   \code{\link{split}(x, index)} should be a valid expression.
#' @param f
#'   A \link{function} or \link{list} of one or more functions to be applied
#'   to subsets of \code{x}. These are recycled to match the number of levels
#'   of \code{index} (unused levels are not dropped). Each function must
#'   accept an initial \link{vector} argument matching \code{\link{typeof}(x)}
#'   (or \code{typeof(x[[j]])} for all \code{j} if \code{x} is a
#'   \link[=data.frame]{data frame}) and return a vector of the same length
#'   (but not necessarily of the same type). It is expected that returned
#'   vectors have compatible \link{mode}s: if one has an \link{atomic} mode,
#'   then all should have an atomic mode (likewise for modes \code{"list"}
#'   and \code{"expression"}).
#'
#' @details
#' Let \code{k = split(seq_len(f), f)}.
#' For \link{vector}s \code{x}, function \code{f[[i]]} is applied
#' to \code{x[k[[i]]]}.
#' For \link[=data.frame]{data frame}s \code{x}, function \code{f[[i]]}
#' is applied to \code{x[[j]][k[[i]]]} for all \code{j}.
#'
#' @return
#' If \code{x} is a \link{vector}, then a vector of the same length.
#' If \code{x} is a \link[=data.frame]{data frame}, then a data frame
#' of the same length, with the same number of rows.
#'
#' @examples
#' ## x <- 1:10
#' ## lpapply(x, index = gl(2L, 5L), f = list(cumprod, function(x) x - mean(x)))
#' ##
#' ## x <- as.data.frame(replicate(3L, c(exp(rnorm(10L)), qlogis(runif(10L)))))
#' ## lpapply(x, index = gl(2L, 10L), f = list(log, plogis))
#'
#' @keywords internal
lpapply <- function(x, index, f) {
  if (!is.list(f)) {
    f <- list(f)
  }
  if (is.data.frame(x)) {
    do_call <- function(g, y) { y[] <- lapply(y, g); y }
  } else {
    do_call <- function(g, y) { g(y) }
  }
  split(x, index) <- Map(do_call, y = split(x, index), g = f)
  x
}

#' Compute Wald confidence intervals
#'
#' Computes approximate confidence intervals around a maximum likelihood
#' estimate (assuming that the estimator is approximately normally distributed).
#'
#' @param estimate,se
#'   \link[=numeric]{Numeric} vectors of equal length supplying maximum
#'   likelihood estimates of parameters and corresponding approximate
#'   standard errors.
#' @param level
#'   A number in the interval (0,1) indicating a confidence level.
#'
#' @details
#' Confidence limits are computed as
#' \code{estimate[i] + c(-1, 1) * sqrt(q) * se[i]},
#' where \code{q = \link{qchisq}(level, df = 1)}.
#' See the \href{https://en.wikipedia.org/wiki/Wald_test}{Wald test}.
#'
#' @return
#' A \link{numeric} \link{matrix} with \code{length(estimate)} rows
#' and 2 columns \code{lower} and \code{upper} giving approximate
#' confidence limits.
#'
#' @examples
#' ## estimate <- rep_len(0, 6L)
#' ## se <- exp(rnorm(6L, sd = 0.1))
#' ## do_wald(estimate = estimate, se = se, level = 0.95)
#'
#' @keywords internal
#' @importFrom stats qchisq
do_wald <- function(estimate, se, level) {
  q <- qchisq(level, df = 1)
  n <- length(estimate)
  lu <- estimate + rep.int(sqrt(q) * c(-1, 1), c(n, n)) * se
  dim(lu) <- c(n, 2L)
  colnames(lu) <- c("lower", "upper")
  lu
}

