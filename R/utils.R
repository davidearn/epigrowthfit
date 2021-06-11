#' Enumerate duplicated strings
#'
#' Replaces the \code{i}th instance of string \code{s}
#' in a \link{character} vector with \code{\link{sprintf}("\%s[\%d]", s, i)}.
#'
#' @param x
#'   For \code{enum_dupl_string}, a \link{character} vector.
#'   For \code{enum_dupl_names}, an \R object with a \code{\link{names}}
#'   \link[=attributes]{attribute}.
#'
#' @return
#' \code{enum_dupl_string(x)} returns an enumerated copy of \code{x}.
#' \code{enum_dupl_names(x)} returns a copy of \code{x} with
#' \code{names(x)} replaced by \code{enum_dupl_string(names(x))}.
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
#'   A \link{logical} scalar. If \code{TRUE}, then \code{\link{NA}} are ignored.
#' @param tol
#'   A positive number. \link[=double]{Double} vectors \code{x} are considered
#'   constant if and only if the distance from \code{min(x)} to \code{max(x)}
#'   is less than \code{tol}.
#'
#' @details
#' \code{TRUE} is returned if \code{length(x) = 0}.
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
#' The result of \code{\link{sweep}(sweep(cor, 1L, sd, `*`), 2L, sd, `*`)},
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

#' Recursively merge lists
#'
#' Recursively merge elements of a partially specified \link{list}
#' into a completely specified list of defaults.
#'
#' @param x
#'   A \link{list}.
#' @param template
#'   A \link{list} into which `x` is merged.
#'
#' @details
#' Recursion proceeds as follows:
#'
#' If \code{x} is \code{\link{NULL}}, then \code{template} is replaced
#' with \code{\link{NULL}}.
#'
#' Otherwise, if \code{x} and \code{template} are both \link{list}s
#' or are both \link{atomic} vectors and \code{template} does not have
#' \link{names}, then \code{x} is recycled to the length of \code{template}.
#'
#' Otherwise, if \code{x} and \code{template} are both \link{list}s and
#' \code{template} does not have \link{names}, then elements of \code{x}
#' are recursively merged into the so-named elements of \code{template}.
#'
#' Otherwise, if \code{x} and \code{template} are both \link{atomic} vectors
#' and \code{template} has \link{names}, then elements of \code{x} replace
#' the so-named elements of \code{template}.
#'
#' Otherwise, \code{x} is discarded and \code{template} is preserved as is.
#'
#' @return
#' A \link{list}.
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
#' # rlmerge(x, template)
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
#' \link[=split]{Split}s a supplied \link[=data.frame]{data frame}
#' on one or more variables in turn.
#'
#' @param x
#'   A \link[=data.frame]{data frame}.
#' @param by
#'   A subset of \code{\link{seq_along}(x)} or \code{\link{names}(x)}
#'   indicating (in order) variables on which to split \code{x}.
#' @param drop
#'   A \link{logical} scalar passed to \code{\link{split}}.
#'
#' @return
#' A recursive \link{list} of \link[=data.frame]{data frame}s
#' with \code{length(by)} "layers", unless \code{length(by) = 0},
#' in which case \code{x} is returned as is.
#'
#' @examples
#' f <- function() sample(letters[1:5], 20L, replace = TRUE)
#' d <- as.data.frame(replicate(3L, f()))
#' # rsplit(d, 1:3)
#'
#' @keywords internal
rsplit <- function(x, by = integer(0L), drop = FALSE) {
  if (length(by) == 0L) {
    return(x)
  }
  l <- split(x, x[[by[1L]]], drop = drop)
  lapply(l, rsplit, by = by[-1L], drop = drop)
}
