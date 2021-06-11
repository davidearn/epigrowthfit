#' Better axis
#'
#' A replacement for [graphics::axis()] allowing the user
#' to define the extent of the axis line.
#'
#' @param a,b
#'   Lower and upper limits of the axis line in user coordinates.
#'   The default (`NULL`) is to use values from `par("usr")`.
#' @param side,at,labels,...
#'   Arguments to [graphics::axis()].
#'
#' @return
#' A numeric vector. The value of `axis(side, at, labels, ...)`.
#'
#' @importFrom graphics axis
baxis <- function(side, a = NULL, b = NULL,
                  at = NULL, labels = TRUE, ...) {
  if (is.null(a) || is.null(b)) {
    gp <- par("usr", "xlog", "ylog")
    if (trunc(side[1L]) %% 2 == 1) {
      usr <- gp$usr[1:2]
      if (gp$xlog) {
        usr <- 10^usr
      }
    } else {
      usr <- gp$usr[3:4]
      if (gp$ylog) {
        usr <- 10^usr
      }
    }
    if (is.null(a)) {
      a <- usr[1L]
    }
    if (is.null(b)) {
      b <- usr[2L]
    }
  }
  dots <- list(...)
  bak <- dots$lwd.ticks
  dots$lwd.ticks <- 0
  l <- list(side = side, at = c(a, b), labels = c("", ""))
  do.call(axis, c(l, dots))

  dots$lwd <- 0
  dots$lwd.ticks <- bak
  l <- list(side = side, at = at, labels = labels)
  do.call(axis, c(l, dots))
}

#' Date axis
#'
#' Adds an axis to the bottom of the current plot and labels
#' it with day, month, and year, taking care to ensure that
#' labels are nicely spaced.
#'
#' @param origin
#'   A Date scalar. It is assumed that horizontal user coordinates
#'   measure time as a number of days since time 00:00:00 on date
#'   `origin`.
#' @param plot
#'   A logical scalar. If `FALSE`, then an axis is not drawn.
#'   Useful when the return value is desired but graphical
#'   output is not.
#' @param minor,major
#'   Named lists of arguments to [graphics::axis()], affecting the
#'   appearance of the minor (day or month) and major (month or year)
#'   axes, respectively.
#'
#' @return
#' A list of numeric vectors `minor` and `major` giving
#' the positions of minor and major axis labels in user
#' coordinates.
#'
#' @details
#' The date axis consists of minor and major axes. The content
#' of these axes depends entirely on `w = diff(par("usr"))[1:2]`,
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
#' @keywords internal
#' @importFrom graphics axis par
daxis <- function(origin = .Date(0), plot = TRUE,
                  minor = NULL, major = NULL) {
  usr <- par("usr")
  t0 <- min(ceiling(usr[1L]), floor(usr[2L]))
  t1 <- max(ceiling(usr[1L]), floor(usr[2L]), t0 + 1)
  d0 <- origin + t0
  d1 <- origin + t1
  w <- t1 - t0

  ## Determine tick coordinates and labels
  if (w <= 210) {
    ## Days
    by <- c(1, 2, 4, 7, 14)[w <= c(14, 28, 56, 112, 210)][1L]
    at_minor_as_date <- seq(d0, d1, by = by)
    at_minor <- julian(at_minor_as_date, origin = d0)
    labels_minor <- ymd(at_minor_as_date, 3L)

    ## Months
    if (ymd(d0, 2L) == ymd(d1, 2L)) {
      at_major_as_date <- d0
      at_major <- 0
    } else {
      at_major_as_date <- seq(dceiling(d0, "month"), d1, by = "month")
      at_major <- julian(at_major_as_date, origin = d0)
      if (at_major[1L] > w / 8) {
        at_major_as_date <- c(d0, at_major_as_date)
        at_major <- c(0, at_major)
      }
    }
    labels_major <- months(at_major_as_date, abbreviate = TRUE)
    any_major <- TRUE
  }
  else if (w <= 3 * 365) {
    ## Months
    by <- c(1L, 2L, 3L)[w <= c(1, 2, 3) * 365][1L]
    at_minor_as_date <- seq(dceiling(d0, "month"), d1, by = paste(by, "months"))
    at_minor <- julian(at_minor_as_date, origin = d0)
    labels_minor <- months(at_minor_as_date, abbreviate = TRUE)

    ## Years
    if (ymd(d0, 1L) == ymd(d1, 1L)) {
      at_major_as_date <- d0
      at_major <- 0
    } else {
      at_major_as_date <- seq(dceiling(d0, "year"), d1, by = "year")
      at_major <- julian(at_major_as_date, origin = d0)
      if (at_major[1L] > w / 8) {
        at_major_as_date <- c(d0, at_major_as_date)
        at_major <- c(0, at_major)
      }
    }
    labels_major <- ymd(at_major_as_date, 1L)
    any_major <- TRUE
  } else {
    ## Years
    by <- ceiling(ceiling(w / 365) / 7)
    at_minor_as_date <- seq(dceiling(d0, "year"), d1 + (by + 1) * 365, by = paste(by, "years"))
    at_minor <- julian(at_minor_as_date, origin = d0)
    labels_minor <- ymd(at_minor_as_date, 1L)
    at_minor <- c(at_minor, (at_minor[-1L] + at_minor[-length(at_minor)]) / 2)
    length(labels_minor) <- length(at_minor)
    any_major <- FALSE
  }

  if (plot) {
    ## Minor axis
    l <- list(
      side = 1,
      at = t0 + at_minor,
      labels = labels_minor
    )
    do.call(baxis, c(l, minor))

    ## Major axis
    if (any_major) {
      l <- list(
        side = 1,
        at = t0 + at_major,
        labels = labels_major
      )
      do.call(baxis, c(l, major))
    }
  }

  list(
    minor = t0 + at_minor,
    major = if (any_major) t0 + at_major else numeric(0L)
  )
}

#' Get nicely formatted tick labels
#'
#' Generates nice "mantissa x 10^power" tick labels for axes,
#' at least for count data.
#'
#' @param at
#'   A numeric vector listing tick positions in user coordinates,
#'   probably generated by [graphics::axTicks()].
#'
#' @return
#' An expression vector of length `length(at)` listing tick labels.
#'
#' @noRd
get_yax_labels <- function(at) {
  ## Exponential notation split into mantissa and power
  mp <- matrix(unlist(strsplit(sprintf("%.6e", at), "e")), ncol = 2L, byrow = TRUE)

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
    labels <- parse(text = sprintf("%s %%*%% 10^%s", man, pow))
  }
  if (0 %in% at) {
    labels[at == 0] <- 0
  }
  labels
}

#' Find size for y-axis tick labels
#'
#' Find `cex.axis` such that the widest _y_-axis tick label generated by
#' `axis(labels, cex.axis, las = 1)` is exactly `mex` margin lines wide.
#'
#' @param labels
#'   A character or expression vector listing _y_-axis tick labels.
#' @param mex
#'   A positive number. The (desired) width of the widest tick label
#'   as a number of margin lines.
#' @param cex,font,csi
#'   Graphical parameters affecting the width of a string as a number
#'   of margin lines. The default (`NULL`) is to obtain values using
#'   [graphics::par()].
#'
#' @return
#' `mex * csi / mlw`, where `mlw` is the width of the widest tick label
#' in inches given `cex` and `font`.
#'
#' @keywords internal
#' @importFrom graphics par strwidth
get_yax_cex <- function(labels, mex, cex = NULL, font = NULL, csi = NULL) {
  gp <- par(c("cex", "font", "csi"))
  if (is.null(cex)) {
    cex <- gp$cex
  }
  if (is.null(font)) {
    font <- gp$font
  }
  if (is.null(csi)) {
    csi <- gp$csi
  }
  mex * csi / max(strwidth(labels, units = "inches", cex = cex, font = font))
}

add_height_to_user <- function(h, y0, log) {
  if (log) y0 * 10^h else y0 + h
}

#' @importFrom graphics par
add_lines_to_user <- function(n, y0, log) {
  gp <- par(c("cxy", "mai", "mar", "pin", "usr"))
  if (log) {
    ## Height of top margin in user coordinates
    u_m3 <- gp$mai[3L] * (gp$usr[4L] - gp$usr[3L]) / gp$pin[2L]
    ## Height of one line in user coordinates
    u_1l <- u_m3 / gp$mar[3L]
    return(y0 * 10^(n * u_1l))
  }
  y0 + n * gp$cxy[2L]
}

#' @importFrom graphics par
user_to_lines <- function(y, log) {
  gp <- par(c("cxy", "mai", "mar", "pin", "usr"))
  if (log) {
    ## Height of top margin in user coordinates
    u_m3 <- gp$mai[3L] * (gp$usr[4L] - gp$usr[3L]) / gp$pin[2L]
    ## Height of one line in user coordinates
    u_1l <- u_m3 / gp$mar[3L]
    return((log10(y) - gp$usr[4L]) / u_1l)
  }
  (y - gp$usr[4L]) / gp$cxy[2L]
}
