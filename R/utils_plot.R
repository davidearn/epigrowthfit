#' Define plot options
#'
#' Sets parameters controlling the graphical output of \code{\link{plot.egf}}.
#'
#' @param points
#'   A named \link{list} of arguments to \code{\link{points}},
#'   affecting the appearance of observed data.
#'   [\code{type != "rt2"} only.]
#' @param points_short,points_long
#'   Alternatives to \code{points} used for counts over intervals
#'   shorter or longer than \code{dt} days.
#'   [\code{type = "interval"} only.]
#' @param lines
#'   A named \link{list} of arguments to \code{\link{lines}},
#'   affecting the appearance of predicted curves.
#'   [\code{type != "rt2"} only.]
#' @param polygon
#'   A named \link{list} of arguments to \code{\link{polygon}},
#'   affecting the appearance of confidence bands on predicted curves.
#'   [\code{type != "rt2"} only.]
#' @param rect
#'   A named \link{list} of arguments to \code{\link{rect}},
#'   affecting the appearance of fitting windows.
#'   [\code{type != "rt2"} only.]
#' @param rect_bg_panel,rect_bg_plab
#'   A named \link{list} of arguments to \code{\link{rect}},
#'   affecting the appearance of panel backgrounds and panel title underlays.
#'   [\code{type = "rt2"} only.]
#' @param abline
#'   A named \link{list} of arguments to \code{\link{abline}},
#'   affecting the appearance of the line drawn at \code{y = 0}.
#'   [\code{type = "rt1"} only.]
#' @param segments
#'   A named \link{list} of arguments to \code{\link{segments}},
#'   affecting the appearance of the line segments drawn at \code{y = r}
#'   \code{x$model$curve} is \code{"exponential"}, \code{"logistic"},
#'   or \code{"richards"}.
#'   [\code{type = "rt1"} only.]
#' @param axis_x_Date_minor,axis_x_Date_major,axis_x_numeric,axis_y
#'   Named \link{list}s of arguments to \code{\link{axis}},
#'   affecting the appearance of plot axes. \code{axis_x_Date_*}
#'   are used for \code{time_as = "Date"}. \code{axis_x_numeric}
#'   is used for \code{time_as = "numeric"}.
#' @param box
#'   A named \link{list} of arguments to \code{\link{box}},
#'   affecting the appearance of the box drawn around the plot region.
#'   [\code{type != "rt2"} only.]
#' @param title_main,title_sub,title_xlab,title_ylab,title_plab
#'   Named \link{list}s of arguments to \code{\link{title}},
#'   affecting the appearance of plot, axis, and panel titles.
#' @param text_tdoubling_estimate,text_tdoubling_ci,text_tdoubling_caption
#'   Named \link{list}s of arguments to \code{\link{text}},
#'   affecting the appearance of initial doubling times printed
#'   in the top margin.
#'   [\code{type != "rt2"} only.]
#' @param colorRamp
#'   A named \link{list} of arguments to \code{\link{colorRamp}},
#'   defining heat map colour palette.
#'
#' @details
#' Unsupported and unmodifiable options are silently discarded.
#' Modifiable options that are unspecified are assigned a default
#' value defined internally. A complete list of defaults can be
#' retrieved with \code{l <- egf_plot_control()}.
#'
#' Setting an argument to \code{\link{NULL}} has the effect of
#' suppressing the corresponding plot element. For example,
#' to suppress observed data, set \code{points} to \code{NULL},
#' and perhaps also \code{points_short} and \code{points_long}.
#'
#' @return
#' A named list.
#'
#' @export
egf_plot_control <- function(points, points_short, points_long,
                             lines, polygon, rect,
                             rect_bg_panel, rect_bg_plab,
                             abline, segments,
                             axis_x_Date_minor, axis_x_Date_major,
                             axis_x_numeric, axis_y, box,
                             title_main, title_sub,
                             title_xlab, title_ylab, title_plab,
                             text_tdoubling_estimate, text_tdoubling_ci,
                             text_tdoubling_caption, colorRamp) {
  nf <- names(formals(egf_plot_control))
  values <- mget(nf, ifnotfound = NA, inherits = FALSE)
  defaults <- list(
    points = list(
      pch = 21,
      col = "#BBBBBB",
      bg = "#DDDDDD",
      cex = 1
    ),
    points_short = list(
      pch = 1,
      col = "#882255",
      bg = NA,
      cex = 1
    ),
    points_long = list(
      pch = 16,
      col = "#882255",
      bg = NA,
      cex = 1
    ),
    lines = list(
      lty = 1,
      lwd = 2.5,
      col = "#44AA99"
    ),
    polygon = list(
      col = "#44AA9960",
      border = NA,
      lty = 1,
      lwd = 1
    ),
    rect = list(
      col = "#DDCC7740",
      border = NA,
      lty = 1,
      lwd = 1
    ),
    rect_bg_panel = list(
      col = "black",
      border = NA,
      lty = 1,
      lwd = 1
    ),
    rect_bg_panel = list(
      col = "#00000080",
      border = NA,
      lty = 1,
      lwd = 1
    ),
    abline <- list(
      lty = 2,
      lwd = 1,
      col = "black"
    ),
    segments <- list(
      lty = 3,
      lwd = 2,
      col = "black"
    ),
    axis_x_Date_minor = list(
      mgp = c(3, 0.25, 0),
      lwd = 1,
      col = "black",
      lwd.ticks = 1,
      col.ticks = "black",
      tcl = -0.2,
      gap.axis = 0,
      col.axis = "black",
      cex.axis = 0.9,
      font.axis = 1,
      family = "sans",
      xpd = FALSE
    ),
    axis_x_Date_major = list(
      mgp = c(3, 1.25, 0),
      lwd = 1,
      col = "black",
      lwd.ticks = 1,
      col.ticks = "black",
      tcl = 0,
      gap.axis = 0,
      col.axis = "black",
      cex.axis = 1.2,
      font.axis = 1,
      family = "sans",
      xpd = FALSE
    ),
    axis_x_numeric = list(
      mgp = c(3, 0.7, 0),
      lwd = 1,
      col = "black",
      lwd.ticks = 1,
      col.ticks = "black",
      tcl = -0.5,
      gap.axis = 0,
      col.axis = "black",
      cex.axis = 0.9,
      font.axis = 1,
      family = "sans",
      xpd = FALSE
    ),
    axis_y = list(
      mgp = c(3, 0.7, 0),
      lwd = 1,
      col = "black",
      lwd.ticks = 1,
      col.ticks = "black",
      tcl = -0.5,
      gap.axis = NA,
      col.axis = "black",
      cex.axis = 0.9,
      font.axis = 1,
      family = "sans",
      xpd = FALSE
    ),
    box = list(
      bty = "l",
      lty = 1,
      lwd = 1,
      col = "black"
    ),
    title_main = list(
      adj = 0,
      col.main = "black",
      cex.main = 1,
      font.main = 2,
      family = "sans"
    ),
    title_sub = list(
      adj = 0,
      col.sub = "black",
      cex.sub = 0.75,
      font.sub = 2,
      family = "sans"
    ),
    title_xlab = list(
      adj = 0.5,
      col.lab = "black",
      cex.lab = 1,
      font.lab = 1,
      family = "sans"
    ),
    title_ylab = list(
      adj = 0.5,
      col.lab = "black",
      cex.lab = 1,
      font.lab = 1,
      family = "sans"
    ),
    title_plab <- list(
      col.lab = "white",
      cex.lab = 1,
      font.lab = 2,
      family = "sans"
    ),
    text_tdoubling_estimate = list(
      col = "black",
      cex = 0.7,
      font = 2,
      family = "sans"
    ),
    test_tdoubling_ci = list(
      col = "black",
      cex = 0.7,
      font = 1,
      family = "sans"
    ),
    text_tdoubling_caption = list(
      col = "black",
      cex = 0.7,
      font = 1,
      family = "sans"
    ),
    colorRamp <- list(
      colors = c(
        "#364B9A", "#4A7BB7", "#6EA6CD", "#98CAE1",
        "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366",
        "#F67E4B", "#DD3D2D", "#A50026"
      ),
      bias = 1,
      space = "rgb",
      interpolate = "linear"
    )
  )

  lmerge <- function(value, default) {
    if (is.null(value)) {
      return(NULL)
    }
    if (!is.list(value) || length(value) == 0L) {
      return(default)
    }
    m <- match(names(value), names(default), 0L)
    default[m] <- value[m > 0L]
    default
  }
  control <- Map(lmerge, value = values, default = defaults)
  class(control) <- c("egf_plot_control", "list")
  control
}


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
#' A \link{numeric} vector.
#' This is the value of \code{\link{axis}(side, at, labels, ...)}.
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
#'   A \link{Date} scalar. It is assumed that horizontal user coordinates
#'   measure time as a number of days since time 00:00:00 on Date \code{origin}.
#' @param do_plot
#'   A \link{logical} scalar. If \code{FALSE}, then an axis is not drawn.
#'   Useful when the return value is desired but graphical output is not.
#' @param minor,major
#'   Named \link{list}s of arguments to \code{\link{axis}}, affecting
#'   the appearance of the minor (day or month) and major (month or year)
#'   axes, respectively.
#'
#' @return
#' A \link{list} of \link{numeric} vectors \code{minor} and \code{major}
#' giving the positions of minor and major axis labels in user coordinates.
#'
#' @details
#' The date axis consists of minor and major axes.
#' The content of these axes depends entirely on
#' \code{w = \link{diff}(\link{par}("usr"))[1:2]},
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
    labels_minor <- ymd(at_minor_as_date, "d")

    ## Months
    if (ymd(d0, "m") == ymd(d1, "m")) {
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
    if (ymd(d0, "y") == ymd(d1, "y")) {
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
    labels_major <- ymd(at_major_as_date, "y")
    any_major <- TRUE
  } else {
    ## Years
    by <- ceiling(ceiling(w / 365) / 7)
    at_minor_as_date <- seq(dceiling(d0, "year"), d1 + (by + 1) * 365, by = paste(by, "years"))
    at_minor <- julian(at_minor_as_date, origin = d0)
    labels_minor <- ymd(at_minor_as_date, "y")
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
#' Generates nice \code{"mantissa x 10^power"} tick labels for axes,
#' at least for count data.
#'
#' @param at
#'   A \link{numeric} vector listing tick positions in user coordinates,
#'   probably generated by \code{\link{axis}}.
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
#' Find \code{cex.axis} such that the widest \eqn{y}-axis tick label
#' generated by \code{\link{axis}(labels, cex.axis, las = 1)} is exactly
#' \code{mex} margin lines wide.
#'
#' @param labels
#'   A \link{character} or \link{expression} vector listing \eqn{y}-axis
#'   tick labels.
#' @param mex
#'   A positive number. The (desired) width of the widest tick label
#'   as a number of margin lines.
#' @param cex,font,csi
#'   Graphical parameters affecting the width of a string as a number
#'   of margin lines. The default (\code{\link{NULL}}) is to obtain values
#'   using \code{\link{par}}.
#'
#' @return
#' \code{mex * csi / mlw}, where \code{mlw} is the width
#' of the widest tick label in inches given \code{cex} and \code{font}.
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
