#' Better axis
#'
#' A replacement for \code{\link{axis}} allowing the user to define
#' the extent of the axis line.
#'
#' @param side,at,labels,...
#'   Arguments to \code{\link{axis}}.
#' @param a,b
#'   Lower and upper limits of the axis line in user coordinates.
#'   The default (\code{NULL}) is to take the limit from
#'   \code{\link{par}("usr")}.
#'
#' @return
#' A numeric vector identical to that returned by
#' \code{\link{axis}(side, at, labels, \dots)}.
#'
#' @noRd
#' @importFrom graphics axis par
baxis <- function(side, a = NULL, b = NULL, at = NULL, labels = TRUE, ...) {
    stopifnot((side <- as.integer(side)) %in% 1:4)
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
#' Adds an axis to the current plot and labels it with day, month, and year,
#' taking care to ensure that labels are nicely spaced.
#'
#' @param side
#'   An integer indicating a side of the plot on which to draw the axis,
#'   passed to \code{\link{axis}}.
#' @param origin
#'   A \link{Date}. Horizontal user coordinates measure time
#'   as a number of days since \code{\link{julian}(origin)} days
#'   after 1970-01-01 00:00:00.
#' @param minor,major
#'   Named \link{list}s of arguments to \code{\link{axis}}, affecting the
#'   appearance of the minor (day or month) and major (month or year) axes,
#'   respectively. Passing \code{NULL} suppresses the corresponding axis.
#'
#' @return
#' A list of numeric vectors \code{minor} and \code{major} giving the positions
#' of minor and major axis labels in user coordinates.
#'
#' @details
#' The \link{Date} axis is a calendar consisting of minor and major axes.
#' The content of these axes depends entirely on
#' \code{w = diff(\link{par}("usr"))[1:2]},
#' the difference between the extreme times in days.
#'
#' If \code{w <= 210} (7 months), then days are placed on the minor
#' axis, months are placed on the major axis, and years are not
#' shown. The spacing on the minor axis in days is the value of
#' \code{c(1, 2, 4, 7, 14)[w <= c(14, 28, 56, 112, 210)][1]}.
#' Hence, for example, if \code{w} is greater than 112 and less than
#' or equal to 210, then the spacing is 14-daily.
#'
#' Otherwise, if \code{w <= 3*365} (3 years), then months are placed
#' on the minor axis, years are placed on the major axis, and days
#' are not shown. The spacing on the minor axis in months is the
#' value of \code{c(1, 2, 3)[w <= c(1, 2, 3) * 365][1]}.
#'
#' Otherwise, if \code{w > 3*365} (3 years), then years are placed on
#' the minor axis, and days and months are not shown. The spacing on the
#' minor axis is the value of \code{ceiling(ceiling(w / 365) / 7)}.
#'
#' @noRd
#' @importFrom graphics axis par
Daxis <- function(side, origin = .Date(0), minor = list(), major = list()) {
    stopifnot((side <- as.integer(side)) %in% 1:4,
              inherits(origin, "Date"),
              length(origin) == 1L,
              is.finite(origin),
              is.list(minor) || is.null(minor),
              is.list(major) || is.null(major))
    usr <- par("usr")[if (side %% 2L == 1L) 1:2 else 3:4]
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

        major_at <- double(0L)
        major <- NULL
    }

    ## Minor axis
    if (!is.null(minor)) {
        args <- list(side = side, at = t0 + minor_at, labels = minor_labels)
        do.call(baxis, c(args, minor))
    }
    ## Major axis
    if (!is.null(major)) {
        args <- list(side = side, at = t0 + major_at, labels = major_labels)
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
#' An \link{expression} vector of length \code{length(at)} listing tick labels.
#'
#' @noRd
get_scientific_labels <- function(at) {
    ## Exponential notation split into mantissa and power
    mp <- matrix(unlist1(strsplit(sprintf("%.6e", at), "e")),
                 ncol = 2L, byrow = TRUE)

    ## Greatest number of digits after mantissa decimal,
    ## ignoring trailing zeros
    digits <- max(nchar(sub("0+$", "", mp[, 1L]))) - 2L

    ## Mantissa reformatted with exactly 'digits' digits after decimal
    man <- sprintf("%.*e", digits, as.double(mp[, 1L]))

    ## Power reformatted without leading zeros
    pow <- as.character(as.double(mp[, 2L]))

    ## Format nonzero labels as "mantissa x 10^power".
    ## Shorten to "10^power" if nonzero mantissas are all 1.
    ## Use "0" if mantissa is 0.
    if (all(as.double(man) %in% c(0, 1))) {
        labels <- parse(text = sprintf("10^%s", pow))
    } else {
        labels <- parse(text = sprintf("%s %%*%% 10^%s", man, pow))
    }
    if (0 %in% at) {
        labels[at == 0] <- expression(0)
    }
    labels
}

#' Calculate a space-filling text size
#'
#' Finds the character expansion factor (what is multiplied by
#' \code{\link{par}("cex")} to obtain the actual magnification)
#' necessary for text to span a given width or height.
#'
#' @param text
#'   A character or \link{expression} vector.
#'   If the length is greater than 1, then only the string or expression
#'   with the maximum width or height is used to calculate the result.
#' @param target
#'   A positive number indicating a target width or height.
#' @param units
#'   A \link{character} string indicating the units of \code{target}.
#' @param horizontal
#'   A \link{logical} flag. If \code{TRUE}, \code{target} represents
#'   a target width rather than a target height.
#' @param ...
#'   Graphical parameters passed
#'   to \code{\link{strwidth}} or \code{\link{strheight}}.
#'
#' @return
#' A positive number.
#'
#' @noRd
#' @importFrom graphics strwidth strheight grconvertX grconvertY
get_sfcex <- function(text, target, units = c("lines", "inches", "user"),
                      horizontal = TRUE, ...) {
    if (length(text) == 0L) {
        return(1)
    }
    measure <- if (horizontal) strwidth else strheight
    inches_current <- max(measure(text, units = "inches", ...))
    if (inches_current == 0) {
        return(1)
    }
    convert <- if (horizontal) grconvertX else grconvertY
    inches_target <- target * diff(convert(c(0, 1), match.arg(units), "inches"))
    inches_target / inches_current
}

#' Modify colour transparency
#'
#' A drop-in replacement for \code{scales::alpha}.
#'
#' @param colour
#'   A numeric or character vector listing colours; see \code{\link{col2rgb}}.
#' @param alpha
#'   A numeric vector with elements in the interval [0,1] listing alpha channel
#'   values.
#'
#' @return
#' A character vector listing colours with indicated transparency.
#'
#' @noRd
#' @importFrom grDevices col2rgb rgb
alpha <- function(colour, alpha) {
    m <- t(col2rgb(colour, alpha = FALSE))
    rgb(m[, 1:3, drop = FALSE], alpha = 255 * alpha, maxColorValue = 255)
}
