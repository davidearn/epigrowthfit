##' Decompose a Date-Time Object
##'
##' Extracts year, month, and day (in UTC) from a date-time object.
##'
##' @param x
##'   a \code{Date}, \code{POSIXct}, or \code{POSIXlt} object.
##' @param which
##'   a character string composed of the letters \samp{y m d},
##'   indicating which of (and in what order) year, month, and day
##'   should be returned.
##' @param drop
##'   a logical.  If \code{drop = TRUE}, then the matrix result
##'   is coerced to a vector when one of its dimensions is less than 2.
##'
##' @return
##' \code{X[, j, drop = drop]}, where \code{X} is an integer matrix
##' with \code{length(x)} rows and 3 columns listing year, month, and
##' day, and \code{j} is the index vector specified by \code{which}.
##'
##' @examples
##' x <- .Date(sample.int(10000L, 10L))
##' ymd(x)
##' ymd(x, which = "d")

ymd <-
function(x, which = "ymd", drop = TRUE) {
	stopifnot(inherits(x, c("Date", "POSIXct", "POSIXlt")),
	          is.character(which),
	          length(which) == 1L,
	          !is.na(which),
	          is.logical(drop),
	          length(drop) == 1L,
	          !is.na(drop))
	x <- as.POSIXct(x)
	attr(x, "tzone") <- "UTC"
	x <- as.POSIXlt(x)
	ans <- c(1900L + x$year, 1L + x$mon, x$mday)
	dim(ans) <- c(length(x), 3L)
	dimnames(ans) <- list(names(x), c("y", "m", "d"))
	j <- match(strsplit(which, "")[[1L]], c("y", "m", "d"), 0L)
	ans[, j, drop = drop]
}

##' Round a Date Vector
##'
##' Rounds each element of a \link{Date} vector to the previous or next
##' midnight (\code{00:00:00}),
##' midnight on a first-of-the-month (\code{YYYY-MM-01}), or
##' midnight on a first-of-the-year (\code{YYYY-01-01}).
##'
##' @param x a \code{Date} vector.
##' @param to a character string.
##'
##' @return
##' \code{x} with elements rounded.
##'
##' @examples
##' n <- 10L
##' x <- .Date(sample.int(10000L, n) + rnorm(n))
##' dmy <- c("day", "month", "year")
##'
##' l <- sapply(dmy, .ceiling.Date, x = x, simplify = FALSE)
##' ceiling.x <- data.frame(unclass.x = unclass(x), x = x, l)
##'
##' l <- sapply(dmy, .floor.Date, x = x, simplify = FALSE)
##' floor.x <- replace(ceiling.x, dmy, l)

.ceiling.Date <-
function(x, units = c("day", "month", "year")) {
	stopifnot(inherits(x, "Date"))
	units <- match.arg(units)
	if (units == "day")
		return(.Date(ceiling(unclass(x))))
	x <- as.POSIXlt(x)
	off <- x$mday > 1L | x$hour > 0L | x$min > 0L | x$sec > 0
	if (units == "month") {
		x$mon <- x$mon + off
		x$year <- x$year + (off <- x$mon > 11L)
		x$mon[off] <- 0L
	}
	else {
		x$year <- x$year + (x$mon > 0L | off)
		x$mon[] <- 0L
	}
	x$mday[] <- 1L
	as.Date(x)
}

.floor.Date <-
function(x, units = c("day", "month", "year")) {
	stopifnot(inherits(x, "Date"))
	units <- match.arg(units)
	if (units == "day")
		return(.Date(floor(unclass(x))))
	x <- as.POSIXlt(x)
	x$mday[] <- 1L
	if (units == "year")
		x$mon[] <- 0L
	as.Date(x)
}
