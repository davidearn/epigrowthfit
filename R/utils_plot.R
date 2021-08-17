#' Better axis
#'
#' A replacement for \code{\link{axis}} allowing the user to define
#' the extent of the axis line.
#'
#' @param a,b
#'   Lower and upper limits of the axis line in user coordinates.
#'   The default (\code{\link{NULL}}) is to take the limit from
#'   \code{\link{par}("usr")}.
#' @param side,at,labels,...
#'   Arguments to \code{\link{axis}}.
#'
#' @return
#' A \link{numeric} vector identical to
#' \code{\link{axis}(side, at, labels, \dots)}.
#'
#' @importFrom graphics axis
baxis <- function(side, a = NULL, b = NULL, at = NULL, labels = TRUE, ...) {
  side <- 1L + (as.integer(side) - 1L) %% 4L
  dots <- list(...)

  if (is.null(a) || is.null(b)) {
    gp <- par("usr", "xlog", "ylog")
    if (side %% 2L == 1L) {
      usr <- gp$usr[1:2]
      log <- gp$xlog
    } else {
      usr <- gp$usr[3:4]
      log <- gp$ylog
    }
    if (log) {
      usr <- 10^usr
    }
    if (is.null(a)) {
      a <- usr[1L]
    }
    if (is.null(b)) {
      b <- usr[2L]
    }
  }

  args <- list(side = side, at = c(a, b), labels = c("", ""))
  do.call(axis, c(args, replace(dots, "lwd.ticks", list(0))))

  args <- list(side = side, at = at, labels = labels)
  do.call(axis, c(args, replace(dots, c("lwd", "lwd.ticks"), list(0, dots[["lwd.ticks"]]))))
}

#' Date axis
#'
#' Adds an axis to the bottom of the current plot and labels
#' it with day, month, and year, taking care to ensure that
#' labels are nicely spaced.
#'
#' @param origin
#'   A \link{Date}. It is assumed that horizontal user coordinates
#'   measure time as a number of days since the time
#'   \code{\link{unclass}(origin)} days after 1970-01-01 00:00:00.
#' @param minor,major
#'   Named \link{list}s of arguments to \code{\link{axis}}, affecting
#'   the appearance of the minor (day or month) and major (month or year)
#'   axes, respectively.
#' @param show_minor,show_major
#'   \link[=logical]{Logical} flags.
#'   If \code{FALSE}, then the corresponding axis is not drawn.
#'
#' @return
#' A \link{list} of \link{numeric} vectors \code{minor} and \code{major}
#' giving the positions of minor and major axis labels in user coordinates.
#'
#' @details
#' The \link{Date} axis is a calendar consisting of minor
#' and major axes. The content of these axes depends entirely
#' on \code{w = \link{diff}(\link{par}("usr"))[1:2]},
#' the difference between the extreme times in days.
#'
#' If \code{w <= 210} (7 months), then days are placed on the minor
#' axis, months are placed on the major axis, and years are not
#' shown. The spacing on the minor axis in days is the value of
#' \code{c(1, 2, 4, 7, 14)[w <= c(14, 28, 56, 112, 210)][1]}.
#' Hence, for example, if `w` is greater than 112 and less than
#' or equal to 210, then the spacing is 14-daily.
#'
#' Otherwise, if \code{w <= 3*365} (3 years), then months are placed
#' on the minor axis, years are placed on the major axis, and days
#' are not shown. The spacing on the minor axis in months is the
#' value of \code{c(1, 2, 3)[w <= c(1, 2, 3) * 365][1]}.
#'
#' Otherwise, if \code{w > 3*365} (3 years), then years are placed on
#' the minor axis, and days and months are not shown. The spacing on the
#' minor axis is the value of \code{\link{ceiling}(ceiling(w / 365) / 7)}.
#'
#' @keywords internal
#' @importFrom graphics axis par
Daxis <- function(origin = .Date(0), minor = NULL, major = NULL,
                  show_minor = TRUE, show_major = TRUE) {
  usr <- par("usr")[1:2]
  Dusr <- origin + usr
  D0 <- min(Dceiling(Dusr[1L]), Dfloor(Dusr[2L]))
  D1 <- max(Dceiling(Dusr[1L]), Dfloor(Dusr[2L]), D0 + 1)
  t0 <- julian(D0, origin = origin)
  t1 <- julian(D1, origin = origin)
  w <- t1 - t0

  ## Determine tick coordinates and labels
  if (w <= 210) {
    ## Days
    by <- c(1, 2, 4, 7, 14)[w <= c(14, 28, 56, 112, 210)][1L]
    minor_at_as_Date <- seq(D0, D1, by = by)
    minor_at <- julian(minor_at_as_Date, origin = D0)
    minor_labels <- ymd(minor_at_as_Date, "d")

    ## Months
    if (ymd(D0, "m") == ymd(D1, "m")) {
      major_at_as_Date <- D0
      major_at <- 0
    } else {
      major_at_as_Date <- seq(Dceiling(D0, "m"), D1, by = "m")
      major_at <- julian(major_at_as_Date, origin = D0)
      if (major_at[1L] > w / 8) {
        major_at_as_Date <- c(D0, major_at_as_Date)
        major_at <- c(0, major_at)
      }
    }
    major_labels <- months(major_at_as_Date, abbreviate = TRUE)
  }
  else if (w <= 3 * 365) {
    ## Months
    by <- c(1L, 2L, 3L)[w <= c(1, 2, 3) * 365][1L]
    minor_at_as_Date <- seq(Dceiling(D0, "m"), D1, by = paste(by, "m"))
    minor_at <- julian(minor_at_as_Date, origin = D0)
    minor_labels <- months(minor_at_as_Date, abbreviate = TRUE)

    ## Years
    if (ymd(D0, "y") == ymd(D1, "y")) {
      major_at_as_Date <- D0
      major_at <- 0
    } else {
      major_at_as_Date <- seq(Dceiling(D0, "y"), D1, by = "y")
      major_at <- julian(major_at_as_Date, origin = D0)
      if (major_at[1L] > w / 8) {
        major_at_as_Date <- c(D0, major_at_as_Date)
        major_at <- c(0, major_at)
      }
    }
    major_labels <- ymd(major_at_as_Date, "y")
  } else {
    ## Years
    by <- ceiling(ceiling(w / 365) / 7)
    minor_at_as_Date <- seq(Dceiling(D0, "year"), D1 + (by + 1) * 365, by = paste(by, "y"))
    minor_at <- julian(minor_at_as_Date, origin = D0)
    minor_labels <- ymd(minor_at_as_Date, "y")
    minor_at <- c(minor_at, (minor_at[-1L] + minor_at[-length(minor_at)]) / 2)
    length(minor_labels) <- length(minor_at)

    major_at <- numeric(0L)
    show_major <- FALSE
  }

  ## Minor axis
  if (show_minor) {
    args <- list(
      side = 1,
      at = t0 + minor_at,
      labels = minor_labels
    )
    do.call(baxis, c(args, minor))
  }
  ## Major axis
  if (show_major) {
    args <- list(
      side = 1,
      at = t0 + major_at,
      labels = major_labels
    )
    do.call(baxis, c(args, major))
  }
  list(minor = t0 + minor_at, major = t0 + major_at)
}

#' Get nicely formatted tick labels
#'
#' Generates nice \code{"mantissa x 10^power"} tick labels for axes,
#' at least for count data.
#'
#' @param at
#'   A \link{numeric} vector listing tick positions in user coordinates,
#'   probably generated by \code{\link{axTicks}}.
#'
#' @return
#' An \link{expression} vector of length \code{\link{length}(at)}
#' listing tick labels.
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
  ## Shorten to "10^power" if nonzero mantissas are all 1.
  ## Use "0" if mantissa is 0.
  if (all(as.numeric(man) %in% c(0, 1))) {
    labels <- parse(text = sprintf("10^%s", pow))
  } else {
    labels <- parse(text = sprintf("%s %%*%% 10^%s", man, pow))
  }
  if (0 %in% at) {
    labels[at == 0] <- expression(0)
  }
  labels
}

#' Find size for y-axis tick labels
#'
#' Find \code{cex.axis} such that the widest \eqn{y}-axis tick label
#' generated by \code{\link{axis}(labels, cex.axis, las = 1, \dots)}
#' spans a desired number of margin lines.
#'
#' @param labels
#'   A \link{character} or \link{expression} vector listing \eqn{y}-axis
#'   tick labels.
#' @param mex
#'   A positive number. The desired maximum label width as a number
#'   of margin lines.
#' @param ...
#'   Graphical parameters passed to \code{\link{strwidth}}.
#'
#' @return
#' A positive number.
#'
#' @keywords internal
#' @importFrom graphics par strwidth
get_yax_cex <- function(labels, mex, ...) {
  mex * par("csi") / max(strwidth(labels, units = "inches", ...))
}
