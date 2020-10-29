#' Annotate a date axis
#'
#' @description
#' Labels day, month, and year on a horizontal time axis, taking care
#' to ensure that labels are nicely spaced. Intended mainly for internal
#' use within [plot.egf_init()] and [`plot.egf()`].
#'
#' @param left Left endpoint of the axis in user coordinates.
#' @param right Right endpoint of the axis in user coordinates.
#' @param refdate A Date scalar. `t0` and `t1` represent a number of
#'   days since this date.
#' @param tcl A numeric scalar. Passed to [graphics::axis()] argument
#'   `tcl` when creating the minor axis. Defines tick length.
#'   See [graphics::par()].
#' @param line A numeric vector of length 2. Elements are passed
#'   to [graphics::axis()] argument `mgp` (the second component)
#'   when creating the minor and major axes, respectively. Defines
#'   the distance between the axis and the top of tick labels.
#'   See [`par()`][graphics::par()].
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
#' The date axis consists of minor and major axes. The content
#' of these axes depends greatly on `w = round(right-left)`,
#' roughly the difference between the first and last date in days.
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
#' All nonzero lengths of `line`, `cex.axis`, and `col.axis`
#' are tolerated. Vectors of length 1 are recycled. Only the
#' first two elements of vectors of length 3 or greater are used.
#'
#' @keywords internal
#' @export
#' @importFrom graphics axis
daxis <- function(left, right, refdate,
                  tcl = -0.2,
                  line = c(0.05, 1),
                  col.axis = c("black", "black"),
                  cex.axis = c(0.7, 0.85)) {
  if (!is.numeric(tcl) || length(tcl) != 1) {
    stop("`tcl` must be numeric and have length 1.")
  }
  if (!is.numeric(line) || length(line) == 0) {
    stop("`line` must be numeric and have nonzero length.")
  }
  if (is.numeric(col.axis) || is.character(col.axis)) {
    if (length(col.axis) == 0) {
      stop("`col.axis` must have nonzero length.")
    } else if (is.numeric(col.axis) && isTRUE(any(col < 0))) {
      stop("Elements of `col.axis` must be non-negative.")
    }
  } else {
    stop("`col.axis` must be a numeric or character vector.")
  }
  if (!is.numeric(cex.axis) || length(cex.axis) == 0) {
    stop("`cex.axis` must be numeric and have nonzero length.")
  }
  for (a in c("line", "col.axis", "cex.axis")) {
    val <- get(a)
    if (length(val) == 1) {
      assign(a, rep(val, 2))
    }
  }


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
  minor_axis <- function() {
    axis(side = 1, at = t0 + at, labels = labels,
         tcl = tcl, mgp = c(3, line[1], 0),
         cex.axis = cex.axis[1], col.axis = col.axis[1])
    invisible(NULL)
  }
  major_axis <- function() {
    axis(side = 1, at = t0 + at, labels = labels,
         tick = FALSE, mgp = c(3, line[2], 0),
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
