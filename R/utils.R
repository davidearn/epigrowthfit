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
#' @name enum_dupl_string
#' @noRd
NULL

enum_dupl_string <- function(x) {
  f <- factor(x)
  n <- tabulate(f)
  i <- unsplit(lapply(n, seq_len), f)
  sprintf("%s[%d]", x, i)
}

enum_dupl_names <- function(x) {
  if (!is.null(names(x))) {
    names(x) <- enum_dupl_string(names(x))
  }
  x
}

#' Correlation to covariance matrix conversion
#'
#' Perform the inverse of \code{\link{cov2cor}}.
#'
#' @param cor
#'   A \link{numeric} \link{matrix}, typically (but not necessarily)
#'   symmetric positive definite.
#' @param sd
#'   A \link{numeric} vector of length \code{\link{nrow}(cor)}.
#'
#' @return
#' The result of
#' \code{\link{sweep}(\link{sweep}(cor, 1L, sd, `*`), 2L, sd, `*`)},
#' obtained slightly more efficiently, following \code{\link{cov2cor}}.
#'
#' @examples
#' X <- replicate(6L, rnorm(10L))
#' V <- cov(X, X)
#' all.equal(V, cor2cov(cov2cor(V), sqrt(diag(V))))
#'
#' @noRd
cor2cov <- function(cor, sd) {
  stopifnot(
    is.matrix(cor),
    is.numeric(cor),
    is.numeric(sd),
    length(sd) == nrow(cor)
  )
  cor[] <- sd * cor * rep(sd, each = length(sd))
  cor
}

#' Concatenate, wrap, collapse
#'
#' Concatenates \R objects, formats the resulting string in paragraphs
#' using \code{\link{strwrap}}, and collapses paragraph lines with newlines.
#' Under default settings, the result is a string that prints in the console
#' as one or more nicely wrapped paragraphs.
#'
#' @param ...
#'   Zero or more \R objects (typically \link{atomic} vectors of length 1),
#'   to be coerced to \link{character} and concatenated with no separator.
#' @param width
#'   A positive integer passed to \code{\link{strwrap}},
#'   indicating a target column width.
#'
#' @return
#' A \link{character} vector of length 1.
#'
#' @noRd
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
#' \code{function(x) inverse.rle(literal_rle(x))}
#' is an identity function for atomic \code{x}, and
#' \code{function(x) literal_rle(inverse.rle(x))}
#' is an identity function for lists \code{x} generated
#' by \code{literal_rle}.
#'
#' @return
#' A \link{list} similar to \code{\link{rle}(x)},
#' though not inheriting from \link{class} \code{"rle"}.
#'
#' @examples
#' x <- rep.int(c(0, NA, NaN, 1), 1:4)
#' rle_x <- literal_rle(x)
#' identical(x, inverse.rle(rle_x))
#'
#' @noRd
literal_rle <- function(x) {
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
#' x <- c(NA, NA, 1, NA, 2, 2, 3, NA)
#' locf(x)
#' locf(x, x0 = 0)
#'
#' @noRd
locf <- function(x, x0 = NULL) {
  if (!anyNA(x)) {
    return(x)
  }
  rle_x <- literal_rle(x)
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

#' Apply length-preserving functions to ragged vectors
#'
#' Modifies ragged vectors in place by replacing each group of elements
#' with the result of applying a length-preserving function to that group
#' of elements.
#'
#' @param x
#'   A \link{vector} or \link[=data.frame]{data frame}.
#' @param index
#'   A \link{factor} (insofar that \code{\link{as.factor}(index)} is a factor)
#'   defining a grouping of the elements or rows of \code{x}.
#' @param f
#'   A \link{function} or \link{list} of one or more functions to be applied
#'   to the subsets of \code{x} defined by \code{index}. Functions are recycled
#'   to the number of levels of \code{index} (unused levels are not dropped).
#'   Each function must accept an initial \link{vector} argument matching
#'   \code{\link{typeof}(x)} (or \code{typeof(x[[j]])} for all \code{j}
#'   if \code{x} is a \link[=data.frame]{data frame}) and return a vector
#'   of the same length (though not necessarily of the same type). It is
#'   expected that returned vectors have compatible \link{mode}s so that
#'   they can be concatenated in a sensible way.
#'
#' @details
#' Let \code{f} be a list of \code{\link{nlevels}(index)} functions,
#' and let \code{k = \link{split}(\link{seq_along}(index), index)}.
#' For vectors \code{x}, function \code{f[[i]]} is applied to
#' \code{x[k[[i]]]}.
#' For data frames \code{x}, function \code{f[[i]]} is applied to
#' \code{x[[j]][k[[i]]]} for all \code{j}.
#'
#' @return
#' If \code{x} is a \link{vector}, then a vector of the same length.
#' If \code{x} is a \link[=data.frame]{data frame}, then a data frame
#' with the same dimensions.
#'
#' @examples
#' x <- 1:10
#' in_place_ragged_apply(x, index = gl(2L, 5L), f = list(cumprod, function(x) x - mean(x)))
#'
#' x <- as.data.frame(replicate(3L, c(exp(rnorm(5L)), qlogis(runif(5L)))))
#' in_place_ragged_apply(x, index = gl(2L, 5L), f = list(log, plogis))
#'
#' @noRd
in_place_ragged_apply <- function(x, index, f) {
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
#' A \link{numeric} \link{matrix} with \code{\link{length}(estimate)} rows
#' and 2 columns \code{lower} and \code{upper} giving approximate confidence
#' limits.
#'
#' @examples
#' estimate <- rnorm(6L, 0, 1)
#' se <- rlnorm(6L, 0, 0.1)
#' do_wald(estimate = estimate, se = se, level = 0.95)
#'
#' @noRd
#' @importFrom stats qchisq
do_wald <- function(estimate, se, level) {
  q <- qchisq(level, df = 1)
  n <- length(estimate)
  lu <- estimate + rep.int(sqrt(q) * c(-1, 1), c(n, n)) * se
  dim(lu) <- c(n, 2L)
  colnames(lu) <- c("lower", "upper")
  lu
}

