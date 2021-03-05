#' Get number of days since a date
#'
#' Subtracts a Date scalar from a Date vector and returns the
#' result as an integer vector instead of a "difftime" object.
#'
#' @param date A Date vector.
#' @param since A Date scalar.
#'
#' @return
#' An integer vector giving the number of days between `date`
#' and `since`. Equal to `as.integer(date - since)`.
#'
#' @noRd
days <- function(date, since = as.Date("1970-01-01")) {
  as.integer(date - since)
}

#' Difference a vector of dates
#'
#' Takes differences of a Date vector and returns the
#' result as an integer vector instead of a "difftime" object.
#'
#' @param date A Date vector.
#' @param lag,differences Arguments to [diff.Date()].
#'
#' @return
#' An integer vector listing the desired differences in days.
#' Equal to `as.integer(diff(date, lag, differences))`.
#'
#' @noRd
ddiff <- function(date, lag = 1L, differences = 1L) {
  as.integer(diff(date, lag, differences))
}

#' Decompose a vector of dates
#'
#' Extracts year, month, and day from a Date vector.
#'
#' @param date
#'   A Date vector.
#' @param which
#'   A subset of `1:3` indicating which of year, month, and day
#'   should be returned.
#' @param drop
#'   A logical scalar. If `drop = TRUE` and one of `date` and `which`
#'   has length less than 2, then an integer vector is returned instead
#'   of a matrix.
#'
#' @return
#' `X[, which, drop]`, where `X` is an integer matrix with
#' `length(date)` rows and 3 columns listing year, month, and day.
#'
#' @noRd
ymd <- function(date, which = 1:3, drop = TRUE) {
  if (length(date) == 0L || length(which) == 0L) {
    x <- integer(0L)
    if (!drop) {
      dim(x) <- c(length(date), length(which))
    }
    return(x)
  }
  X <- matrix(as.integer(unlist(strsplit(as.character(date), "-"))),
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
#' @param date A Date vector.
#' @param to A character string.
#'
#' @return
#' A Date vector of length `length(date)`.
#' If `to = "month"`, then elements are firsts-of-the-month (`"YYYY-MM-01"`).
#' If `to = "year"`, then elements are firsts-of-the-year (`"YYYY-01-01"`).
#'
#' @noRd
dceiling <- function(date, to = c("month", "year")) {
  if (length(date) == 0L) {
    x <- integer(0L)
    class(x) <- "Date"
    return(x)
  }

  to <- match.arg(to)
  X <- as.data.frame(ymd(date, drop = FALSE))

  if (to == "month") {
    X$m <- X$m + (X$d > 1L)
    X$y <- X$y + (X$m == 13L)
    X$m[X$m == 13L] <- 1L
    as.Date(paste(X$y, X$m, "1", sep = "-"))
  } else { # year
    X$y <- X$y + (X$m > 1L || X$d > 1L)
    as.Date(paste(X$y, "1", "1", sep = "-"))
  }
}

#' Annotate a date axis
#'
#' Labels day, month, and year on the bottom axis of a plot,
#' taking care to ensure that labels are nicely spaced.
#'
#' @param left
#'   Left endpoint of the bottom axis in user coordinates.
#' @param right
#'   Right endpoint of the bottom axis in user coordinates.
#' @param refdate
#'   A Date scalar. `left` and `right` represent numbers
#'   of days since this date.
#' @param plot
#'   A logical scalar. If `FALSE`, then the axis is not drawn.
#'   Useful if only the return value is desired.
#' @param tcl,mgp2,col.axis,cex.axis,font.axis
#'   Arguments to [graphics::axis()], which are recycled to length 2.
#'   (`mgp2` is passed as the second component of argument `mgp`.)
#'   The first and second elements control the appearance of the minor
#'   and major axes, respectively. Documentation can be found under
#'   [graphics::par()].
#'
#' @return
#' An integer vector listing positions of minor axis labels.
#'
#' @details
#' `daxis()` assumes that horizontal user coordinates measure
#' time in days since `refdate`.
#'
#' The date axis consists of minor and major axes. The content
#' of these axes depends entirely on `w = floor(right)-ceiling(left)`,
#' the difference between the extreme times in days.
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
#' Otherwise, if `w > 3*365` (3 years), then years are placed on
#' the minor axis, and days and months are not shown. The spacing
#' of the minor axis is the value of `ceiling(ceiling(w / 365) / 7)`.
#'
#' @noRd
#' @importFrom graphics axis par
daxis <- function(left = NULL, right = NULL, refdate, plot = TRUE,
                  tcl = NULL, mgp2 = NULL,
                  col.axis = NULL, cex.axis = NULL, font.axis = NULL) {
  if (is.null(left)) {
    left <- par("usr")[1L]
  }
  if (is.null(right)) {
    right <- par("usr")[2L]
  }
  t0 <- min(ceiling(left), floor(right))
  t1 <- max(ceiling(left), floor(right), t0 + 1)
  d0 <- refdate + t0
  d1 <- refdate + t1
  w <- t1 - t0

  ## Determine tick coordinates and labels
  if (w <= 210) {
    ## Days
    by <- c(1L, 2L, 4L, 7L, 14L)[w <= c(14, 28, 56, 112, 210)][1L]
    at_minor_as_date <- seq(d0, d1, by = by)
    at_minor <- days(at_minor_as_date, since = d0)
    labels_minor <- ymd(at_minor_as_date, 3L)

    ## Months
    if (ymd(d0, 2L) == ymd(d1, 2L)) {
      at_major_as_date <- d0
      at <- 0L
    } else {
      at_major_as_date <- seq(dceiling(d0, "month"), d1, by = "month")
      at_major <- days(at_major_as_date, since = d0)
      if (at_major[1L] > w / 8) {
        at_major_as_date <- c(d0, at_major_as_date)
        at_major <- c(0L, at_major)
      }
    }
    labels_major <- months(at_major_as_date, abbreviate = TRUE)
    any_major <- TRUE
  }
  else if (w <= 3 * 365) {
    ## Months
    by <- c(1L, 2L, 3L)[w <= c(1, 2, 3) * 365][1L]
    at_minor_as_date <- seq(dceiling(d0, "month"), d1, by = paste(by, "months"))
    at_minor <- days(at_minor_as_date, since = d0)
    labels_minor <- months(at_minor_as_date, abbreviate = TRUE)

    ## Years
    if (ymd(d0, 1L) == ymd(d1, 1L)) {
      at_major_as_date <- d0
      at_major <- 0
    } else {
      at_major_as_date <- seq(dceiling(d0, "year"), d1, by = "year")
      at_major <- days(at_major_as_date, since = d0)
      if (at_major[1L] > w / 8) {
        at_major_as_date <- c(d0, at_major_as_date)
        at_major <- c(0L, at_major)
      }
    }
    labels_major <- ymd(at_major_as_date, 1L)
    any_major <- TRUE
  } else {
    ## Years
    nyear <- ceiling(w / 365)
    by <- ceiling(nyear / 7)
    at_minor_as_date <- seq(dceiling(d0, "year"), d1 + (by + 1) * 365, by = paste(by, "years"))
    at_minor <- days(at_minor_as_date, since = d0)
    labels_minor <- ymd(at_minor_as_date, 1L)
    at_minor <- c(at_minor, (at_minor[-1L] + at_minor[-length(at_minor)]) / 2)
    length(labels_minor) <- length(at_minor)
    any_major <- FALSE
  }

  if (plot) {
    if (is.null(tcl)) {
      tcl <- par("tcl")
    }
    tcl <- tcl[1L]
    if (is.null(mgp2)) {
      mgp2 <- par("mgp2")[2L] + c(0, 1)
    }
    mgp2 <- rep_len(mgp2, 2L)
    if (is.null(col.axis)) {
      col.axis <- par("col.axis")
    }
    col.axis <- rep_len(col.axis, 2L)
    if (is.null(cex.axis)) {
      cex.axis <- par("cex.axis")
    }
    cex.axis <- rep_len(cex.axis, 2L)
    if (is.null(font.axis)) {
      font.axis <- par("font.axis")
    }
    font.axis <- rep_len(font.axis, 2L)

    ## Axis line
    axis(
      side = 1L,
      at = c(left, right),
      labels = c("", ""),
      lwd.ticks = 0
    )
    ## Minor axis
    axis(
      side = 1L,
      at = t0 + at_minor,
      labels = labels_minor,
      xpd = TRUE,
      gap.axis = 0,
      lwd = 0,
      lwd.ticks = 1,
      tcl = tcl[1L],
      mgp = c(3, mgp2[1L], 0),
      col.axis = col.axis[1L],
      cex.axis = cex.axis[1L],
      font.axis = font.axis[1L]
    )
    ## Major axis
    if (any_major) {
      axis(
        side = 1L,
        at = t0 + at_major,
        labels = labels_major,
        xpd = TRUE,
        gap.axis = 0,
        lwd = 0,
        mgp = c(3, mgp2[2L], 0),
        col.axis = col.axis[2L],
        cex.axis = cex.axis[2L],
        font.axis = font.axis[2L]
      )
    }
  }

  at_minor
}

#' Get nicely formatted tick labels
#'
#' Generates nice "mantissa x 10^power" tick labels for axes,
#' at least for count data.
#'
#' @param at
#'   A numeric vector listing tick coordinates,
#'   probably generated by [graphics::axTicks()].
#'
#' @return
#' An expression vector of length `length(at)` listing tick labels.
#'
#' @noRd
get_labels <- function(at) {
  ## Exponential notation split into mantissa and power
  mp <- matrix(unlist(strsplit(sprintf("%.6e", at), "e")),
    ncol = 2L,
    byrow = TRUE
  )

  ## Greatest number of digits after mantissa decimal,
  ## ignoring trailing zeros
  digits <- max(nchar(sub("0+$", "", mp[, 1L]))) - 2L

  ## Mantissa reformatted with exactly `digits` digits after decimal
  man <- sprintf("%.*e", digits, as.numeric(mp[, 1L]))

  ## Power reformatted without leading zeros
  pow <- as.character(as.numeric(mp[, 2L]))

  ## Format nonzero labels as "mantissa x 10^power".
  ## Shorten to "10^power" if nonzero mantissa is only ever 1.
  ## Use "0" if mantissa is zero.
  if (all(as.numeric(man) %in% c(0, 1))) {
    labels <- parse(text = sprintf("10^%s", pow))
  } else {
    labels <- parse(text = sprintf("%s %*% 10^%s", man, pow))
  }
  if (0 %in% at) {
    labels[at == 0] <- expression(0)
  }
  labels
}

#' Find size for y-axis tick labels
#'
#' Find `cex.axis` such that the widest y-axis tick label generated
#' by `axis(labels = x, las = 1, cex.axis)` is exactly `mex` margin
#' lines wide.
#'
#' @param x
#'   A character or expression vector listing tick labels.
#' @param mex
#'   A positive number. The width of the widest tick label
#'   as a number of margin lines. See [graphics::par()].
#' @param csi
#'   A positive number. The width of one margin line in inches.
#'   The default is `par("csi")`. See [graphics::par()].
#'
#' @return
#' `mex * csi / mw`, where `mw` is the width of the widest tick label
#' in inches, given by `max(strwidth(x, units = "inches", cex = 1))`.
#'
#' @noRd
#' @importFrom graphics par strwidth
get_cex_axis <- function(x, mex, csi = NULL) {
  if (is.null(csi)) {
    csi <- par("csi")
  }
  mex * csi / max(strwidth(x, units = "inches", cex = 1))
}

add_height_to_user <- function(h, y0, log) {
  if (log) y0 * 10^h else y0 + h
}

#' @importFrom graphics par
add_lines_to_user <- function(n, y0, log) {
  if (log) {
    ## Height of top margin in user coordinates
    u_m3 <- par("mai")[3L] * diff(par("usr")[3:4]) / par("pin")[2L]
    ## Height of one line in user coordinates
    u_1l <- u_m3 / par("mar")[3L]
    y0 * 10^(n * u_1l)
  } else {
    y0 + n * par("cxy")[2L]
  }
}

#' @importFrom graphics par
user_to_lines <- function(y, log) {
  if (log) {
    ## Height of top margin in user coordinates
    u_m3 <- par("mai")[3L] * diff(par("usr")[3:4]) / par("pin")[2L]
    ## Height of one line in user coordinates
    u_1l <- u_m3 / par("mar")[3L]
    (log10(y) - par("usr")[4L]) / u_1l
  } else {
    (y - par("usr")[4L]) / par("cxy")[2L]
  }
}

inv_log10 <- function(x, log) {
  if (log) 10^x else x
}

make_endpoints <- function(object) {
  refdate <- min(object$frame[[1L]])
  time_split <- split(days(object$frame[[1L]], since = refdate), object$index)
  endpoints <- data.frame(
    .t1 = vapply(time_split, min, 0L),
    .t2 = vapply(time_split, max, 0L)
  )
  row.names(endpoints) <- NULL
  attr(endpoints, "refdate") <- refdate
  endpoints
}
