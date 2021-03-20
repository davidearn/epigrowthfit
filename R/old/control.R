#' Get control default
#'
#' Retrieves the default value of the `control` argument of functions
#' in the \pkg{epigrowthfit} package. Currently, the only function with
#' a `control` argument is [plot.egf()].
#'
#' @param f A character string giving the name of a function.
#' @param ... Optional arguments to the function.
#'
#' @return
#' A named list.
#'
#' @export
get_control_default <- function(f, ...) {
  dots <- list(...)
  if (f == "plot.egf") {
    type <- match.arg(dots$type, formals(f)$type)
    if (type == "rt2") {
      l <- list(
        colorRamp = list(colors = c("#364B9A", "#4A7BB7", "#6EA6CD", "#98CAE1",
                                    "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366",
                                    "#F67E4B", "#DD3D2D", "#A50026"),
                         bias = 1, space = "rgb", interpolate = "linear"),
        ips = 0.25
      )
    } else {
      l <- list(
        box = list(bty = "l", lty = 1, lwd = 1, col = 1),
        xax_numeric = list(mgp = c(3, 0.7, 0), tcl = -0.5, gap.axis = 0,
                           col = 1, lwd = 1, col.ticks = 1, lwd.ticks = 1,
                           col.axis = 1, cex.axis = 0.85, font.axis = 1, xpd = TRUE),
        xax_Date_minor = list(mgp = c(3, 0.25, 0), tcl = -0.2, gap.axis = 0,
                              col = 1, lwd = 1, col.ticks = 1, lwd.ticks = 1,
                              col.axis = 1, cex.axis = 0.85, font.axis = 1, xpd = TRUE),
        xax_Date_major = list(mgp = c(3, 1.25, 0), tcl = 0, gap.axis = 0,
                              col = 1, lwd = 0, col.ticks = 1, lwd.ticks = 1,
                              col.axis = 1, cex.axis = 0.85, font.axis = 1, xpd = TRUE),
        yax = list(mgp = c(3, 0.7, 0), tcl = -0.5, gap.axis = NA,
                   col = 1, lwd = 1, col.ticks = 1, lwd.ticks = 1,
                   col.axis = 1, cex.axis = 0.85, font.axis = 1,
                   xpd = TRUE),
        main = list(adj = 0, col.main = 1, cex.main = 1, font.main = 2),
        sub = list(adj = 0, col.main = 1, cex.main = 0.75, font.main = 2),
        xlab = list(adj = 0.5, col.lab = 1, cex.lab = 1, font.lab = 1),
        ylab = list(adj = 0.5, col.lab = 1, cex.lab = 1, font.lab = 1),
        rect = list(col = "#DDCC7740", border = NA, lty = 1, lwd = 1),
        polygon = list(col = "#44AA9960", border = NA, lty = 1, lwd = 1),
        lines = list(lty = 1, lwd = 2.5, col = "#44AA99"),
        points = list(pch = 21, col = "#BBBBBB", bg = "#DDDDDD", cex = 1),
        text_tdoubling_caption = list(col = 1, cex = 0.7, font = 1),
        text_tdoubling_estimate = list(col = 1, cex = 0.7, font = 2),
        text_tdoubling_ci = list(col = 1, cex = 0.7, font = 1)
      )
      if (type == "interval") {
        l <- c(l, list(
          tol = 0,
          points_short = list(pch = 1, col = "#882255", bg = NA, cex = 1),
          points_long = list(pch = 16, col = "#882255", bg = NA, cex = 1),
          text_short_long = list(pos = 3, offset = 0.3, col = "#BBBBBB", cex = 0.7, font = 2)
        ))
      }
      if (type == "rt1") {
        l <- c(l, list(
          abline = list(lty = 2, lwd = 1, col = "black"),
          segments = list(lty = 3, lwd = 2, col = "black")
        ))
      }
    }
  } else {
    stop(sprintf("No control default implemented for function name `%s`.", f))
  }
  l
}

get_clean_control <- function(control, f, ...) {
  default <- get_control_default(f, ...)
  if (is.null(control) || !inherits(control, "list") || is.null(names(control))) {
    return(default)
  }
  for (s in intersect(names(default), names(control))) {
    if (is.null(control[[s]])) {
      default[s] <- list(NULL)
    } else if (inherits(default[[s]], "list") && !is.null(names(default[[s]])) &&
               inherits(control[[s]], "list") && !is.null(names(control[[s]]))) {
      ss <- intersect(names(default[[s]]), names(control[[s]]))
      default[[s]][ss] <- control[[s]][ss]
    } else {
      default[[s]] <- control[[s]]
    }
  }
  default
}
