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
    type <- if (is.null(dots$type)) "" else dots$type
    list(
      box = list(bty = "l", lty = 1, lwd = 1, col = "black"),
      xax = list(tcl = -0.2, mgp2 = c(0.25, 1.25), col.axis = "black", cex.axis = c(0.85, 1.15), font.axis = 1),
      yax = list(tcl = -0.5, mgp2 = 0.7, col.axis = "black", cex.axis = 0.85, font.axis = 1),
      main = list(adj = 0, col.main = "black", cex.main = c(1, 0.75), font.main = 2),
      ylab = list(adj = 0.5, col.lab = "black", cex.lab = 1, font.lab = 1),
      rect = list(col = "#DDCC7740", border = NA, lty = 1, lwd = 1),
      polygon = list(col = "#44AA9960", border = NA, lty = 1, lwd = 1),
      lines = list(lty = 1, lwd = 2.5, col = "#44AA99"),
      points = list(pch = 21, col = "#BBBBBB", bg = "#DDDDDD", cex = 1),
      points_short = list(pch = 1, col = "#882255", bg = NA, cex = 1),
      points_long = list(pch = 16, col = "#882255", bg = NA, cex = 1),
      text_short_long = list(pos = 3, offset = 0.3, col = "#BBBBBB", cex = 0.7, font = 2),
      text_tdoubling = list(col = "black", cex = 0.7, font = c(2, 1, 1)),
      abline = list(lty = 2, lwd = 1, col = "black"),
      segments = list(lty = 3, lwd = 2, col = "black"),
      colorRamp = list(colors = c("#364B9A", "#4A7BB7", "#6EA6CD", "#98CAE1",
                                  "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366",
                                  "#F67E4B", "#DD3D2D", "#A50026"),
                       bias = 1, space = "rgb", interpolate = "linear")
    )
  } else {
    stop(sprintf("No control default implemented for function name `%s`.", f))
  }
}
