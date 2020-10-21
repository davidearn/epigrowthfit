#' Methods for class "egf_sim"
#'
#' @description
#' Methods for "egf_sim" objects returned by [simulate.egf()].
#'
#' @param x An "egf_sim" object.
#' @param inc One of `"cumulative"` and `"interval"`,
#'   indicating a type of incidence to plot.
#' @param alpha A number in \[0,1\] indicating an alpha level.
#'   Determines the transparency of plotted simulations
#'   (0=invisible, 1=opaque). More simulations tend to require
#'   more transparency.
#' @param ... Optional arguments. Used only by the `plot` method
#'   to specify graphical parameters. See [`plot()`][base::plot()]
#'   and [`par()`][graphics::par()]. Currently, only `xlim` and
#'   `ylim` are implemented. Any additional parameters will be
#'   ignored.
#'
#' @return
#' The `plot` method returns `NULL` invisibly.
#'
#' @details
#' ## Plot elements
#'
#' The bottom axis measures the number of days since
#' `x$object$init$date[1]`. The left axis measures
#' interval or cumulative incidence (depending on `inc`).
#'
#' Simulations are obtained as the columns of `x$cum_inc`
#' or `x$int_inc` (depending on `inc`). They are plotted
#' in grey with transparency at the time points given by
#' `x$time` or `x$time[-1]`, respectively.
#'
#' The predicted curve is obtained as `pred$cum_inc`
#' or `pred$int_inc` (depending on `inc`),
#' where `pred = predict(x$object, time = x$time)`.
#' It is plotted in teal without transparency.
#'
#' @seealso [simulate.egf()]
#' @name egf_sim-methods
NULL

#' @rdname egf_sim-methods
#' @export
#' @import graphics
#' @importFrom scales alpha
plot.egf_sim <- function(x, inc = "cumulative", alpha = 0.5, ...) {
  if (!is.character(inc) || length(inc) != 1 ||
        !inc %in% c("interval", "cumulative")) {
    stop("`inc` must be one of \"interval\", \"cumulative\".")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 ||
        !isTRUE(alpha >= 0 && alpha <= 1)) {
    stop("`alpha` must be a number in [0,1].")
  }


  ### SET UP ###########################################################

  ## Optional graphical parameters
  dots <- list(...)

  ## Simulated data
  x$int_inc <- rbind(NA, x$int_inc)

  ## Predicted curves
  pred <- predict(x$object, time = x$time)[c("time", "cum_inc", "int_inc")]
  pred$int_inc <- c(NA, pred$int_inc)

  ## A way to avoid conditional `if (inc = ...) ... else ...`
  varname <- substr(inc, start = 1, stop = 3) # first three characters
  varname <- paste0(varname, "_inc")

  ## Axis titles
  xlab <- paste("days since", as.character(x$object$init$date[1]))
  ylab <- paste(inc, "incidence")

  ## Axis limits (x)
  xmin <- 0
  xmax <- max(x$time, na.rm = TRUE) * 1.04
  xlim <- if ("xlim" %in% names(dots)) dots$xlim else c(xmin, xmax)

  ## Axis limits (y)
  ymin <- 0
  ymax <- max(x[[varname]], pred[[varname]], na.rm = TRUE) * 1.04
  ylim <- if ("ylim" %in% names(dots)) dots$ylim else c(ymin, ymax)

  ## Line styles
  lines_col_sim <- "#BBBBBB"
  lines_col_pred <- "#44AA99"


  ### PLOT #############################################################

  op <- par(mar = c(3, 5, 1, 2), las = 1, mgp = c(3, 0.7, 0))
  on.exit(par(op))
  plot.new()
  plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i")

  ## Simulated data
  for (j in 1:ncol(x[[varname]])) {
    lines(x$time, x[[varname]][, j],
          lwd = 2, col = alpha(lines_col_sim, alpha = alpha))
  }

  ## Predicted curves
  lines(pred$time, pred[[varname]], lwd = 3, col = lines_col_pred)

  ## Axes
  box(bty = "l")
  axis(side = 1)
  axis(side = 2)

  ## Titles
  title(xlab = xlab, line = 2)
  title(ylab = ylab, line = 3.8)

  invisible(NULL)
}
