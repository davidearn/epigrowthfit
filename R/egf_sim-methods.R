#' Methods for class "egf_sim"
#'
#' @description
#' Methods for "egf_sim" objects returned by [simulate.egf()].
#'
#' @param x An "egf_sim" object.
#' @param inc One of `"cumulative"` and `"interval"`,
#'   indicating a type of incidence to plot.
#' @param col_sim,col_pred Character or numeric scalars specifying
#'   colours for simulated incidence time series and predicted incidence
#'   curves, respectively.
#' @param ... Optional arguments. Used only by the `plot` method
#'   to specify graphical parameters. Currently, only `xlim`,
#'   `ylim`, `xlab`, `ylab`, and `main` are implemented. Further
#'   arguments are ignored.
#'   See [`plot()`][graphics::plot()] and [`par()`][graphics::par()]
#'   for a catalogue of graphical parameters.
#'
#' @return
#' The `plot` method returns `NULL` invisibly.
#'
#' @details
#' ## Plot elements
#'
#' The bottom axis measures the number of days since
#' `x$object$init$date[1]`. The left axis measures interval
#' or cumulative incidence (depending on `inc`). The predicted
#' incidence curve is obtained as `pred$cum_inc` or `pred$int_inc`,
#' where `pred = predict(x$object, time = x$time)`. The simulated
#' incidence time series are obtained as the columns of `x$cum_inc`
#' or `x$int_inc`. These are plotted together behind the predicted
#' incidence curve. To help visualize the distribution of simulated
#' incidence at a given time, assign `col_sim` a sufficiently
#' transparent colour.
#'
#' @seealso [simulate.egf()]
#' @name egf_sim-methods
NULL

#' @rdname egf_sim-methods
#' @export
#' @import graphics
plot.egf_sim <- function(x, inc = "cumulative",
                         col_pred = "#44AA99", col_sim = "#BBBBBB66", ...) {
  if (!is.character(inc) || length(inc) != 1 ||
        !inc %in% c("interval", "cumulative")) {
    stop("`inc` must be one of \"interval\", \"cumulative\".")
  }


  ### SET UP ###########################################################

  ## Optional graphical parameters
  dots <- list(...)

  ## Simulated data
  x$int_inc <- rbind(NA, x$int_inc)

  ## Predicted curve
  pred <- predict(x$object, time = x$time)[c("time", "cum_inc", "int_inc")]
  pred$int_inc <- c(NA, pred$int_inc)

  ## A way to avoid conditional `if (inc = ...) ... else ...`
  varname <- substr(inc, start = 1, stop = 3) # first three characters
  varname <- paste0(varname, "_inc")

  ## Titles
  if ("xlab" %in% names(dots)) {
    xlab <- dots$xlab
  } else {
    xlab <- paste("days since", as.character(x$object$init$date[1]))
  }
  if ("ylab" %in% names(dots)) {
    ylab <- dots$ylab
  } else {
    ylab <- paste(inc, "incidence")
  }
  cstr <- x$object$init$curve
  substr(cstr, 1, 1) <- toupper(substr(cstr, 1, 1)) # capitalize first letter
  if ("main" %in% names(dots)) {
    main <- dots$main
  } else {
    main <- paste0(cstr, " model of ", inc, " incidence\n",
                   "(", ncol(x[[varname]]), " simulations)")
  }

  ## Axis limits (x)
  if ("xlim" %in% names(dots)) {
    xlim <- dots$xlim
  } else {
    xmin <- 0
    xmax <- max(x$time) * 1.04
    xlim <- c(xmin, xmax)
  }

  ## Axis limits (y)
  if ("ylim" %in% names(dots)) {
    ylim <- dots$ylim
  } else {
    ymin <- 0
    ymax <- max(x[[varname]], na.rm = TRUE) * 1.04
    ylim <- c(ymin, ymax)
  }


  ### PLOT #############################################################

  op <- par(
    mar = c(4, 5, 2.7, 0.5) + 0.1,
    las = 1,
    mgp = c(3, 0.7, 0)
  )
  plot.new()
  plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i")

  ## Simulated data
  for (j in 1:ncol(x[[varname]])) {
    lines(x$time, x[[varname]][, j], lwd = 2, col = col_sim)
  }

  ## Predicted curves
  lines(pred$time, pred[[varname]], lwd = 3, col = col_pred)

  ## Box
  box(bty = "l")

  ## Axis (x)
  axis(side = 1, cex.axis = 0.85)

  ## Axis (y)
  yax_at <- axTicks(side = 2)
  if (max(yax_at) < 1e05) {
    yax_labels <- TRUE
    digits <- 0
  } else {
    mp <- matrix(unlist(strsplit(sprintf("%.6e", yax_at), "e")),
                 ncol = 2, byrow = TRUE)
    digits <- max(nchar(sub("0+$", "", mp[, 1]))) - 2
    man <- sprintf(paste0("%.", digits, "e"), as.numeric(mp[, 1]))
    pow <- as.character(as.numeric(mp[, 2]))
    if (all(as.numeric(man) %in% c(0, 1))) {
      yax_labels <- parse(text = paste0("10^", pow))
    } else {
      yax_labels <- parse(text = paste0(man, " %*% 10^", pow))
    }
    if (0 %in% yax_at) {
      yax_labels[yax_at == 0] <- expression(0)
    }
  }
  axis(side = 2, at = yax_at, labels = yax_labels,
       cex.axis = if (digits > 1) 0.65 else 0.85)

  ## Titles
  title(xlab = xlab, line = 3)
  title(ylab = ylab, line = 4)
  title(main = main, line = 1, cex.main = 0.9)

  par(op)
  invisible(NULL)
}
