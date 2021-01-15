#' Get control default
#'
#' Retrieves the default value of the `control` argument of functions
#' in the \pkg{epigrowthfit} package. Currently, the only function with
#' a `control` argument is [plot.egf()].
#'
#' @param f A character string giving the name of a function.
#'
#' @return
#' A named list.
#'
#' @export
get_control_default <- function(f) {
  switch(f,
    plot.egf = list(
      box = list(bty = "l", lty = 1, lwd = 1, col = "black"),
      xax = list(tcl = -0.2, mgp2 = c(0.25, 1.25), col.axis = "black", cex.axis = c(0.85, 1.15), font.axis = 1),
      yax = list(tcl = -0.5, mgp2 = 0.7, col.axis = "black", cex.axis = NA, font.axis = 1),
      xlab = list(line = 3, adj = 0.5, col.lab = "black", cex.lab = 1, font.lab = 1),
      ylab = list(line = 4, adj = 0.5, col.lab = "black", cex.lab = 1, font.lab = 1),
      main = list(line = 1, adj = 0.5, col.main = "black", cex.main = 1, font.main = 2),
      points_main = list(pch = 21, col = "#BBBBBB", bg = "#DDDDDD", cex = 1),
      points_short = list(pch = 1, col = "#882255", bg = NA, cex = 1),
      points_long = list(pch = 16, col = "#882255", bg = NA, cex = 1),
      lines = list(lty = 1, lwd = 2.5, col = "#44AA99"),
      window = list(col = "#DDCC7740", border = NA),
      confband = list(col = "#44AA9960", border = NA),
      text_hl = list(pos = 3, offset = 0.3, col = "#BBBBBB", cex = 0.7, font = 2),
      text_dbl = list(x = NA, y = NA, adj = c(0, 0.5), pos = 4, offset = 1, col = "black", cex = 0.7, font = 1)
    ),
    stop("No control default implemented for function name `f`.")
  )
}
