#' Define plot options
#'
#' Sets parameters controlling the graphical output of \code{\link{plot.egf}}.
#' Here, \code{x}, \code{type}, \code{time_as}, and \code{dt} refer to the
#' so-named arguments of \code{\link{plot.egf}}.
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
#'   affecting the appearance of line segments drawn at \code{y = r}
#'   when \code{x$model$curve} is \code{"exponential"}, \code{"logistic"},
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
#' @param mtext_tdoubling_caption,mtext_tdoubling_estimate,mtext_tdoubling_ci
#'   Named \link{list}s of arguments to \code{\link{mtext}},
#'   affecting the appearance of initial doubling times printed
#'   in the top margin.
#'   [\code{type != "rt2"} only.]
#' @param colorRamp
#'   A named \link{list} of arguments to \code{\link{colorRamp}},
#'   defining heat map colour palette.
#'   [\code{type = "rt2"} only.]
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
                             mtext_tdoubling_caption,
                             mtext_tdoubling_estimate, mtext_tdoubling_ci,
                             colorRamp) {
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
    mtext_tdoubling_caption = list(
      col = "black",
      cex = 0.7,
      font = 1,
      family = "sans",
      adj = 1
    ),
    mtext_tdoubling_estimate = list(
      col = "black",
      cex = 0.7,
      font = 2,
      family = "sans"
    ),
    mtext_tdoubling_ci = list(
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
  do.call(axis, c(args, replace(dots, "lwd", list(0))))
}

#' Date axis
#'
#' Adds an axis to the bottom of the current plot and labels
#' it with day, month, and year, taking care to ensure that
#' labels are nicely spaced.
#'
#' @param origin
#'   A \link{Date}. It is assumed that horizontal user coordinates
#'   measure time as a number of days since time 00:00:00 on Date
#'   \code{origin}.
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
  t0 <- min(ceiling(usr[1L]), floor(usr[2L]))
  t1 <- max(ceiling(usr[1L]), floor(usr[2L]), t0 + 1)
  d0 <- origin + t0
  d1 <- origin + t1
  w <- t1 - t0

  ## Determine tick coordinates and labels
  if (w <= 210) {
    ## Days
    by <- c(1, 2, 4, 7, 14)[w <= c(14, 28, 56, 112, 210)][1L]
    minor_at_as_Date <- seq(d0, d1, by = by)
    minor_at <- julian(minor_at_as_Date, origin = d0)
    minor_labels <- ymd(minor_at_as_Date, "d")

    ## Months
    if (ymd(d0, "m") == ymd(d1, "m")) {
      major_at_as_Date <- d0
      major_at <- 0
    } else {
      major_at_as_Date <- seq(dceiling(d0, "month"), d1, by = "month")
      major_at <- julian(major_at_as_Date, origin = d0)
      if (major_at[1L] > w / 8) {
        major_at_as_Date <- c(d0, major_at_as_Date)
        major_at <- c(0, major_at)
      }
    }
    major_labels <- months(major_at_as_Date, abbreviate = TRUE)
  }
  else if (w <= 3 * 365) {
    ## Months
    by <- c(1L, 2L, 3L)[w <= c(1, 2, 3) * 365][1L]
    minor_at_as_Date <- seq(dceiling(d0, "month"), d1, by = paste(by, "months"))
    minor_at <- julian(minor_at_as_Date, origin = d0)
    minor_labels <- months(minor_at_as_Date, abbreviate = TRUE)

    ## Years
    if (ymd(d0, "y") == ymd(d1, "y")) {
      major_at_as_Date <- d0
      major_at <- 0
    } else {
      major_at_as_Date <- seq(dceiling(d0, "year"), d1, by = "year")
      major_at <- julian(major_at_as_Date, origin = d0)
      if (major_at[1L] > w / 8) {
        major_at_as_Date <- c(d0, major_at_as_Date)
        major_at <- c(0, major_at)
      }
    }
    major_labels <- ymd(major_at_as_Date, "y")
  } else {
    ## Years
    by <- ceiling(ceiling(w / 365) / 7)
    minor_at_as_Date <- seq(dceiling(d0, "year"), d1 + (by + 1) * 365, by = paste(by, "years"))
    minor_at <- julian(minor_at_as_Date, origin = d0)
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
#' generated by \code{\link{axis}(labels, cex.axis, las = 1)} spans
#' a desired number of margin lines.
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
