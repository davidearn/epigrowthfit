#' Validate an atomic vector
#'
#' @description
#' Checks the class, length, and value of an atomic vector and
#' optionally throws an error if desired conditions are not met.
#' Simplifies the task of validating atomic vector arguments of
#' functions.
#'
#' @param x An atomic vector.
#' @param ... Zero or more objects coercible to character,
#'   which are pasted without separator to generate
#'   an error or warning message (depending on `action`).
#' @param what A character vector listing object classes.
#'   Check fails if `x` does not inherit from one of these classes.
#' @param not A character vector listing object classes.
#'   Check fails if `x` inherits from one of these classes.
#' @param len A numeric vector. If `length(len) = 1`, then check
#'   fails if `length(x) != len`. If `length(len) > 1`, then check
#'   fails if `length(x) < len[1]` or `length(x) > len[2]`.
#' @param opt An atomic vector listing accepted values for `x[i]`.
#'   Check fails if `!all(x %in% opt)`.
#' @param val A numeric vector indicating endpoints for the interval
#'   of acceptable values for `x[i]`, assuming that `x` is numeric.
#'   Whether the endpoints are included in this interval is determined
#'   by `rel`.
#'   If `length(val) = 1`, then check fails if
#'   `!all(x == val, na.rm = TRUE)`.
#'   If `length(val) > 1` and (for example) `rel = c(">=", "<=")`,
#'   then check fails if
#'   `!all(x >= val[1] & x <= val[2], na.rm = TRUE)`.
#' @param rel A character vector of length 2 such that
#'   `rel[1] %in% c(">", ">=")` and `rel[2] %in% c("<", "<=")`.
#' @param yes A function or list of functions of one argument returning
#'   either `TRUE` or `FALSE`. Check fails if any of these functions
#'   return `FALSE` when applied to `x`.
#' @param no A function or list of functions of one argument returning
#'   either `TRUE` or `FALSE`. Check fails if any of these functions
#'   return `TRUE` when applied to `x`.
#' @param action One of `"stop"`, `"warning"`, and `"nothing"`,
#'   indicating how check failure should be handled.
#'
#' @return
#' `TRUE` if check passes and `FALSE` otherwise.
#'
#' @details
#' Checks are performed in the order of `names(formals(check))`.
#' If `what` fails then `not` won't be checked. If `not` fails,
#' then `len` won't be checked (and so on).
#'
#' Ironically, no checks on the input are performed, so use with
#' care (e.g., if setting `val`, then make sure that `x` is numeric
#' or that `what` contains `"numeric"`).
#'
#' Integer vectors are accepted as numeric vectors when `what`
#' contains `"numeric"`. Hence `check(1:6, what = "numeric")`
#' returns `TRUE`.
#'
#' @keywords internal
#' @export
check <- function(x, ..., what = NULL, not = NULL, len = NULL,
                  opt = NULL, val = NULL, rel = c(">=", "<="),
                  yes = NULL, no = NULL,
                  action = c("stop", "warning", "nothing")) {
  ## Tolerate integer when checking for numeric
  if ("numeric" %in% what) {
    what <- c(what, "integer")
  }
  ## Tolerate function instead of list of functions
  if (is.function(yes)) {
    yes <- list(yes)
  }
  if (is.function(no)) {
    no <- list(no)
  }
  ## Tolerate, e.g., action = "warn"
  action <- match.arg(action)

  passes <-
    ## Check class
    (is.null(what) || inherits(x, what)) &&
    (is.null(not) || !inherits(x, not)) &&
    ## Check length: equals `len` or in `range(len)`
    (is.null(len) ||
       (length(len) == 1 && length(x) == len) ||
       (length(len) > 1 && length(x) >= len[1] && length(x) <= len[2])) &&
    ## Check value: is a subset of `opt`
    (is.null(opt) || all(x %in% opt)) &&
    ## Check value: equals `val` or in `range(val)`
    (is.null(val) ||
       (length(val) == 1 && all(match.fun(rel[1])(x, val[1]), na.rm = TRUE)) ||
       (length(val) > 1 && all(match.fun(rel[2])(x, val[2]), na.rm = TRUE))) &&
    ## Check value: passes functions `yes`
    (is.null(yes) || all(sapply(yes, function(f) f(x)))) &&
    ## Check value: fails functions `no`
    (is.null(no) || all(!sapply(no, function(f) f(x))))
  if (!isTRUE(passes) && action %in% c("stop", "warning")) {
    do.call(action, list(..., call. = FALSE))
  }
  passes
}

#' Annotate a date axis
#'
#' @description
#' Labels day, month, and year on a bottom axis, taking care to ensure
#' that labels are nicely spaced.
#'
#' @param left Left endpoint of the bottom axis in user coordinates.
#' @param right Right endpoint of the bottom axis in user coordinates.
#' @param refdate A Date scalar. `t0` and `t1` represent a number of
#'   days since this date.
#' @param tcl A numeric scalar. Passed to [graphics::axis()] argument
#'   `tcl` when creating the minor axis. Defines tick length.
#'   See [graphics::par()].
#' @param mgp2 A numeric vector of length 2. Elements are passed
#'   to [graphics::axis()] argument `mgp` (the second component)
#'   when creating the minor and major axes, respectively. Defines
#'   the distance between the axis and the top of tick labels.
#'   See [graphics::par()].
#' @param col.axis A numeric or character vector of length 2.
#'   Elements are passed to [graphics::axis()] argument `col.axis`
#'   when creating the minor and major axes, respectively. Defines
#'   the colour of tick labels. See [graphics::par()].
#' @param cex.axis A numeric vector of length 2. Elements are
#'   passed to [graphics::axis()] argument `cex.axis` when creating
#'   the minor and major axes, respectively. Defines the size of
#'   tick labels. See [graphics::par()].
#'
#' @return
#' Returns `NULL` invisibly.
#'
#' @details
#' `daxis()` assumes that horizontal user coordinates measure
#' a number of days since `refdate`.
#'
#' The date axis consists of minor and major axes. The content
#' of these axes depends entirely on `w = round(right-left)`,
#' the rounded difference between the extreme times in days.
#'
#' If `w <= 210` (7 months), then days are placed on the minor
#' axis, months are placed on the major axis, and years are not
#' shown. The spacing of the minor axis in days is the value of
#' `c(1, 2, 4, 7, 14)[w <= c(14, 28, 56, 112, 210)][1]`. Hence,
#' for example, if `w` is greater than 112 and less than or equal
#' to 210, then the spacing is 14-daily.
#'
#' Otherwise, if `w <= 3*365` (3 years), then months are placed
#' on the minor axis, years are placed on the major axis, and days
#' are not shown. The spacing of the minor axis in months is the
#' value of `c(1, 2, 3)[w <= c(1, 2, 3) * 365][1]`.
#'
#' Otherwise (i.e., if `w > 3*365`), years are placed on the minor
#' axis, and days and months are not shown. The spacing of the
#' minor axis is the value of `ceiling(ceiling(w / 365) / 7)`.
#'
#' All nonzero lengths of `mgp2`, `cex.axis`, and `col.axis`
#' are tolerated. Vectors of length 1 are recycled. Only the
#' first two elements of vectors of length 3 or greater are used.
#'
#' @keywords internal
#' @export
#' @importFrom graphics axis
daxis <- function(left, right, refdate,
                  tcl = -0.2,
                  mgp2 = c(0.05, 1),
                  col.axis = "black",
                  cex.axis = c(0.7, 0.85)) {
  ### SET UP ###########################################################

  t0 <- min(ceiling(left), floor(right))
  t1 <- max(ceiling(left), floor(right), t0 + 1)
  d0 <- refdate + t0
  d1 <- refdate + t1
  w <- as.numeric(d1 - d0)
  ymd <- function(date, which = 1:3) {
    cmat <- matrix(unlist(strsplit(as.character(date), "-")),
                   ncol = 3, byrow = TRUE)
    nmat <- apply(cmat, 2, as.numeric)
    if (is.null(dim(nmat))) nmat[which] else nmat[, which]
  }
  ceiling_y <- function(date) {
    y <- ymd(date, 1) + ifelse(ymd(date, 2) == 1 & ymd(date, 3) == 1, 0, 1)
    as.Date(paste(y, 1, 1, sep = "-"))
  }
  ceiling_m <- function(date) {
    y <- ymd(date, 1)
    m <- ymd(date, 2) + ifelse(ymd(date, 3) == 1, 0, 1)
    y <- ifelse(m > 12, y + 1, y)
    m <- ifelse(m > 12, 1, m)
    as.Date(paste(y, m, 1, sep = "-"))
  }
  mgp2 <- rep(mgp2, length.out = 2)
  col.axis <- rep(col.axis, length.out = 2)
  cex.axis <- rep(cex.axis, length.out = 2)
  minor_axis <- function() {
    axis(side = 1, at = t0 + at, labels = labels,
         tcl = tcl, mgp = c(3, mgp2[1], 0), gap.axis = 0,
         cex.axis = cex.axis[1], col.axis = col.axis[1])
    invisible(NULL)
  }
  major_axis <- function() {
    axis(side = 1, at = t0 + at, labels = labels,
         tick = FALSE, mgp = c(3, mgp2[2], 0), gap.axis = 0,
         cex.axis = cex.axis[2], col.axis = col.axis[2])
    invisible(NULL)
  }


  ### PLOT #############################################################

  if (w <= 210) {
    ## Days
    by <- c(1, 2, 4, 7, 14)[w <= c(14, 28, 56, 112, 210)][1]
    at_date <- seq(d0, d1, by = by)
    at <- as.numeric(at_date - d0)
    labels <- ymd(at_date, 3)
    minor_axis()

    ## Months
    if (ymd(d0, 2) == ymd(d1, 2)) {
      at_date <- d0
      at <- 0
    } else {
      at_date <- seq(ceiling_m(d0), d1, by = "month")
      at <- as.numeric(at_date - d0)
      if (at[1] > w / 8) {
        at_date <- c(d0, at_date)
        at <- c(0, at)
      }
    }
    labels <- months(at_date, abbreviate = TRUE)
    major_axis()
  }
  else if (w <= 3 * 365) {
    ## Months
    by <- c(1, 2, 3)[w <= c(1, 2, 3) * 365][1]
    at_date <- seq(ceiling_m(d0), d1, by = paste(by, "months"))
    at <- as.numeric(at_date - d0)
    labels <- months(at_date, abbreviate = TRUE)
    minor_axis()

    ## Years
    if (ymd(d0, 1) == ymd(d1, 1)) {
      at_date <- d0
      at <- 0
    } else {
      at_date <- seq(ceiling_y(d0), d1, by = "year")
      at <- as.numeric(at_date - d0)
      if (at[1] > w / 8) {
        at_date <- c(d0, at_date)
        at <- c(0, at)
      }
    }
    labels <- ymd(at_date, 1)
    major_axis()
  } else {
    ## Years
    nyear <- ceiling(w / 365)
    by <- ceiling(nyear / 7)
    at_date <- seq(ceiling_y(d0), d1 + (by + 1) * 365,
                   by = paste(by, "years"))
    at <- as.numeric(at_date - d0)
    n <- length(at)
    at <- c(at, (at[-1] + at[-n]) / 2)
    labels <- c(ymd(at_date, 1), rep(NA, length(at)-n))
    minor_axis()
  }

  invisible(NULL)
}