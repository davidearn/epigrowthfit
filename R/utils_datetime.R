#' Decompose a Date vector
#'
#' Extracts year, month, and day from a \link{Date} vector.
#'
#' @param x
#'   A \link{Date} vector.
#' @param which
#'   A \link{character} string composed of the letters \samp{y m d},
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
#' # x <- .Date(sample.int(1e4L, 10L))
#' # ymd(x)
#' # ymd(x, which = "d", drop = TRUE)
#'
#' @keywords internal
ymd <- function(x, which = "ymd", drop = TRUE) {
  stop_if_not_string(which)
  stop_if_not_true_false(drop)
  X <- matrix(as.integer(unlist(strsplit(as.character(x), "-"), FALSE, FALSE)),
    nrow = length(x),
    ncol = 3L,
    byrow = TRUE,
    dimnames = list(NULL, c("y", "m", "d"))
  )
  j <- unique(match(strsplit(which, "")[[1L]], colnames(X), 0L))
  X[, j, drop = drop]
}

#' Get ceiling of a date
#'
#' Rounds each Date in a \link{Date} vector to the next first-of-the-month
#' or first-of-the-year.
#'
#' @param x A \link{Date} vector.
#' @param to A \link{character} string.
#'
#' @return
#' \code{x} with elements replaced by firsts-of-the-month (\code{"YYYY-MM-01"})
#' or firsts-of-the-year (\code{"YYYY-01-01"}), depending on \code{to}.
#'
#' @examples
#' # x <- .Date(sample.int(1e4L, 10L))
#' # dceiling(x, to = "month")
#' # dceiling(x, to = "year")
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
