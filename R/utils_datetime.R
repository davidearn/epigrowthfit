#' Decompose a Date vector
#'
#' Extracts year, month, and day from a \link{Date} vector.
#'
#' @param x
#'   A \link{Date} vector.
#' @param which
#'   A character string composed of the letters \samp{y m d},
#'   indicating which of (and in what order) year, month, and day
#'   should be returned.
#' @param drop
#'   A logical scalar. If \code{drop = TRUE}, then the matrix result
#'   is coerced to a vector when one of its dimensions is less than 2.
#'
#' @return
#' \code{X[, j, drop = drop]}, where \code{X} is an \link{integer}
#' \link{matrix} with \code{length(x)} rows and 3 columns listing
#' year, month, and day, and \code{j} is an index vector derived
#' from \code{which}.
#'
#' @examples
#' x <- .Date(sample.int(1e4L, 10L))
#' ymd(x)
#' ymd(x, which = "d", drop = TRUE)
#'
#' @noRd
ymd <- function(x, which = "ymd", drop = TRUE) {
  stopifnot(
    inherits(x, "Date") || is.character(x),
    is.character(which),
    length(which) == 1L,
    !is.na(which),
    is.logical(drop),
    length(drop) == 1L,
    !is.na(drop)
  )
  X <- matrix(NA_integer_, nrow = length(x), ncol = 3L, dimnames = list(names(x), c("y", "m", "d")))
  if (length(x) > 0L) {
    ok <- is.finite(x)
    i <- as.integer(unlist1(strsplit(as.character(x[ok]), "-")))
    X[ok, ] <- matrix(i, nrow = sum(ok), ncol = 3L, byrow = TRUE)
  }
  j <- match(strsplit(which, "")[[1L]], colnames(X), 0L)
  X[, j, drop = drop]
}

#' Round a Date vector
#'
#' Rounds each element of a \link{Date} vector to the previous or next
#' midnight (\code{00:00:00}),
#' midnight on a first-of-the-month (\code{YYYY-MM-01}), or
#' midnight on a first-of-the-year (\code{YYYY-01-01}).
#'
#' @param x A \link{Date} vector.
#' @param to A \link{character} string.
#'
#' @return
#' \code{x} with elements rounded.
#'
#' @examples
#' n <- 10L
#' x <- .Date(sample.int(10000L, n) + rnorm(n))
#' dmy <- c("day", "month", "year")
#'
#' l <- sapply(dmy, Dceiling, x = x, simplify = FALSE)
#' ceiling_x <- data.frame(unclass_x = unclass(x), x = x, l)
#'
#' l <- sapply(dmy, Dfloor, x = x, simplify = FALSE)
#' floor_x <- replace(ceiling_y, dmy, l)
#'
#' @noRd
NULL

Dceiling <- function(x, to = c("day", "month", "year")) {
  stopifnot(inherits(x, "Date"))
  to <- match.arg(to)
  x <- .Date(ceiling(unclass(x)))
  if (to == "day" || !any(ok <- is.finite(x))) {
    return(x)
  }
  X <- as.data.frame(ymd(x[ok], drop = FALSE))
  if (to == "month") {
    X$m <- X$m + (X$d > 1L)
    X$y <- X$y + (i <- X$m == 13L)
    X$m[i] <- 1L
    x[ok] <- as.Date(paste(X$y, X$m, "1", sep = "-"), format = "%Y-%m-%d")
  } else { # to == "year"
    X$y <- X$y + (X$m > 1L || X$d > 1L)
    x[ok] <- as.Date(paste(X$y, "1", "1", sep = "-"), format = "%Y-%m-%d")
  }
  x
}

Dfloor <- function(x, to = c("day", "month", "year")) {
  stopifnot(inherits(x, "Date"))
  to <- match.arg(to)
  x <- .Date(floor(unclass(x)))
  if (to == "day" || !any(ok <- is.finite(x))) {
    return(x)
  }
  X <- as.data.frame(ymd(x[ok], drop = FALSE))
  if (to == "month") {
    x[ok] <- as.Date(paste(X$y, X$m, "1", sep = "-"), format = "%Y-%m-%d")
  } else { # to == "year"
    x[ok] <- as.Date(paste(X$y, "1", "1", sep = "-"), format = "%Y-%m-%d")
  }
  x
}
