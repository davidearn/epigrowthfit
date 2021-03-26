#' Get control default
#'
#' Retrieves the default value of the `control` argument of functions
#' in the \pkg{epigrowthfit} namespace. Currently, the only function
#' with a `control` argument is [plot.egf()].
#'
#' @param f A character string giving the name of a function.
#' @param ... Optional arguments to the function named `f`.
#'
#' @return
#' A named list.
#'
#' @export
get_control_default <- function(f = c("plot.egf"), ...) {
  f <- match.arg(f)
  ff <- formals(f)
  dots <- list(...)
  if (f == "plot.egf") {
    type <- match.arg(dots$type, eval(ff$type))
    time_as <- match.arg(dots$time_as, eval(ff$time_as))
    default <- list(
      box = list(
        bty = "l",
        lty = 1,
        lwd = 1,
        col = 1
      ),
      axis = list(
        x = list(
          Date = list(
            minor = list(
              mgp = c(3, 0.25, 0),
              lwd = 1,
              col = 1,
              lwd.ticks = 1,
              col.ticks = 1,
              tcl = -0.2,
              gap.axis = 0,
              col.axis = 1,
              cex.axis = 0.85,
              font.axis = 1,
              xpd = FALSE
            ),
            major = list(
              mgp = c(3, 1.25, 0),
              lwd = 1,
              col = 1,
              lwd.ticks = 1,
              col.ticks = 1,
              tcl = 0,
              gap.axis = 0,
              col.axis = 1,
              cex.axis = switch(type, rt2 = 1.25, 1.15),
              font.axis = 1,
              xpd = FALSE
            )
          ), # Date
          numeric = list(
            mgp = c(3, 0.7, 0),
            lwd = 1,
            col = 1,
            lwd.ticks = 1,
            col.ticks = 1,
            tcl = -0.5,
            gap.axis = 0,
            col.axis = 1,
            cex.axis = 0.85,
            font.axis = 1,
            xpd = FALSE
          )
        )[[time_as]], # x
        y = list(
          mgp = c(3, 0.7, 0),
          lwd = 1,
          col = 1,
          lwd.ticks = 1,
          col.ticks = 1,
          tcl = -0.5,
          gap.axis = NA,
          col.axis = 1,
          cex.axis = 0.85,
          font.axis = 1,
          xpd = FALSE
        )
      ), # axis
      title = list(
        main = list(
          adj = 0,
          col.main = 1,
          cex.main = 1,
          font.main = 2
        ),
        sub = list(
          adj = 0,
          col.sub = 1,
          cex.sub = 0.75,
          font.sub = 2
        ),
        xlab = list(
          adj = 0.5,
          col.lab = 1,
          cex.lab = 1,
          font.lab = 1
        ),
        ylab = list(
          adj = 0.5,
          col.lab = 1,
          cex.lab = 1,
          font.lab = 1
        )
      ), # title
      rect = list(
        col = "#DDCC7740",
        border = NA,
        lty = 1,
        lwd = 1
      ),
      polygon = list(
        col = "#44AA9960",
        border = NA,
        lty = 1,
        lwd = 1
      ),
      lines = list(
        lty = 1,
        lwd = 2.5,
        col = "#44AA99"
      ),
      points = list(
        pch = 21,
        col = "#BBBBBB",
        bg = "#DDDDDD",
        cex = 1
      ),
      tdoubling = list(
        caption = list(
          col = 1,
          cex = 0.7,
          font = 1
        ),
        estimate = list(
          col = 1,
          cex = 0.7,
          font = 2
        ),
        ci = list(
          col = 1,
          cex = 0.7,
          font = 1
        )
      ),
      special = list(
        tol = 0,
        points = list(
          short = list(
            pch = 1,
            col = "#882255",
            bg = NA,
            cex = 1
          ),
          long = list(
            pch = 16,
            col = "#882255",
            bg = NA,
            cex = 1
          )
        ),
        text = list(
          pos = 3,
          offset = 0.3,
          col = "#BBBBBB",
          cex = 0.7,
          font = 2
        )
      ), # special
      abline <- list(
        lty = 2,
        lwd = 1,
        col = 1
      ),
      segments <- list(
        lty = 3,
        lwd = 2,
        col = "black"
      )
    )

    if (type == "rt2") {
      default <- default[c("axis", "title")]
      default$title$sub <- NULL
      default$title$ylab$adj <- NULL
      default$title$plab <- list(
        col.lab = 0,
        cex.lab = 1,
        font.lab = 2
      )
      default$rect <- list(
        bg = list(
          col = 1,
          border = NA,
          lty = 1,
          lwd = 1
        ),
        plab = list(
          col = "#00000080",
          border = NA,
          lty = 1,
          lwd = 1
        )
      )
      default$colorRamp <- list(
        colors = c(
          "#364B9A", "#4A7BB7", "#6EA6CD", "#98CAE1",
          "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366",
          "#F67E4B", "#DD3D2D", "#A50026"
        ),
        bias = 1,
        space = "rgb",
        interpolate = "linear"
      )
      default$ips <- 0.25
      return(default)
    }
    if (type != "interval") {
      default$special <- NULL
    }
    if (type != "rt1") {
      default$abline <- default$segments <- NULL
    }
    return(default)
  } # f == "plot.egf"
  NULL
}
