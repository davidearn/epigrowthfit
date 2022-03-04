#' Decompose a date-time object
#'
#' Extracts year, month, and day (in UTC) from a date-time object.
#'
#' @param x
#'   A \link{POSIXlt}, \link{POSIXct}, or \link{Date} object.
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
#' x <- .Date(sample.int(1e+04L, 10L))
#' ymd(x)
#' ymd(x, which = "d", drop = TRUE)
#'
#' @noRd
ymd <- function(x, which = "ymd", drop = TRUE) {
    stopifnot(inherits(x, c("Date", "POSIXlt", "POSIXct")),
              is.character(which),
              length(which) == 1L,
              !is.na(which),
              is.logical(drop),
              length(drop) == 1L,
              !is.na(drop))
    x <- as.POSIXct(x)
    attr(x, "tzone") <- "UTC"
    x <- as.POSIXlt(x)
    res <- matrix(c(1900L + x$year, 1L + x$mon, x$mday),
                  nrow = length(x), ncol = 3L,
                  dimnames = list(names(x), c("y", "m", "d")))
    j <- match(strsplit(which, "")[[1L]], colnames(res), 0L)
    res[, j, drop = drop]
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
    if (to == "day") {
        return(.Date(ceiling(unclass(x))))
    }
    x <- as.POSIXlt(x)
    if (to == "month") {
        x$mon <- x$mon + (x$mday > 1L)
        x$year <- x$year + (i <- x$mon > 11L)
        x$mon[i] <- 0L
    } else {
        x$year <- x$year + (x$mon > 0L || x$mday > 1L)
        x$mon[] <- 0L
    }
    x$mday[] <- 1L
    as.Date(x)
}

Dfloor <- function(x, to = c("day", "month", "year")) {
    stopifnot(inherits(x, "Date"))
    to <- match.arg(to)
    if (to == "day") {
        return(.Date(floor(unclass(x))))
    }
    x <- as.POSIXlt(x)
    x$mday[] <- 1L
    if (to == "year") {
        x$mon[] <- 0L
    }
    as.Date(x)
}
