#' Methods for class "egf"
#'
#' @description
#' A method for plotting "egf" objects returned by [egf()].
#'
#' @param x An "egf" object.
#' @param inc One of `"interval"` and `"cumulative"`,
#'   indicating whether to plot `x$cases` (interval incidence)
#'   or `cumsum(x$cases)` (cumulative incidence).
#' @param tol A non-negative number used only if
#'   `inc = "interval"`. `x$ini$cases[i]` is plotted
#'   in blue if `diff(x$init$time)[i] < (1-tol)*m`,
#'   in red if `diff(x$init$time)[i] > (1+tol)*m`,
#'   and in grey otherwise, where `m = median(diff(x$init$time))`.
#'   Assign `Inf` to ensure that everything is grey.
#' @param ... Unused optional arguments.
#'
#' @details
#' If `x$init$time` is not equally spaced, then `inc = "interval"`
#' should be used with caution. `x$init$cases[i]` is the number
#' of cases observed between `x$init$time[i]` and `x$init$time[i+1]`,
#' hence `x$init$cases` roughly scales with `diff(x$init$time)`.
#' Argument `tol` can be used with `inc = "interval"` in order to
#' highlight outliers in `diff(x$init$time)`.
#'
#' @name egf-methods
NULL

#' @rdname egf-methods
#' @export
#' @import graphics
#' @importFrom stats median
plot.egf <- function(x, inc = "interval", tol = 0.025, ...) {
  if (!inherits(x, "egf")) {
    stop("`x` must be an \"egf\" object.")
  }
  if (!is.character(inc) || length(inc) != 1 ||
      !inc %in% c("interval", "cumulative")) {
    stop("`inc` must be an element of `c(\"interval\", \"cumulative\")`.")
  }
  if (inc == "interval") {
    if (!is.numeric(tol) || length(tol) != 1 || !isTRUE(tol >= 0)) {
      stop("`tol` must be a non-negative number.")
    }
  }
  y <- x$init
  dtime <- diff(y$time)
  m <- median(dtime)
  f <- y$first
  l <- y$last
  wgrid <- seq(y$time[f+1], y$time[l+1], by = m)
  op <- par(mar = c(5, 4, 4, 8) + 0.1, las = 1, mgp = c(3, 0.7, 0))
  data <- data.frame(time = y$time[-1], cases = y$cases)
  xlab <- paste("days since", as.character(y$date[1]))
  if (inc == "interval") {
    dtime_min <- (1 - tol) * m
    dtime_max <- (1 + tol) * m
    bg_enum <- 1 + 1 * (dtime < dtime_min) + 2 * (dtime > dtime_max)
    bg <- c("#DDDDDD", "#66CCEE", "#4477AA")[bg_enum]
    plot(cases + 0.1 ~ time, data = data, xlab = xlab,
         log = "y", pch = 21, bg = bg)
    int_inc <- x$eval_int_inc(wgrid)
    lines(wgrid[-1], int_inc, lwd = 2, col = "#EE6677")
  } else if (inc == "cumulative") {
    plot(cumsum(cases) + 0.1 ~ time, data = data, xlab = xlab,
         log = "y", pch = 21, bg = "#DDDDDD")
    cum_inc <- x$eval_cum_inc(wgrid)
    lines(wgrid, cum_inc, lwd = 2, col = "#EE6677")
  }
  abline(v = y$time[c(f,l)+1], lty = 2, col = "#555555")
  axis(side = 3, at = y$time[c(f,l)+1], labels = c(f,l),
       tick = FALSE, mgp = c(3, 0.1, 0))
  mtext("index", side = 3, line = 2)
  ## Model
  mstr <- paste0(y$curve,
                 if (y$include_baseline) "\nbaseline" else "",
                 "\n", y$distr)
  mx <- par("usr")[2] + 0.02 * diff(par("usr")[1:2])
  my <- 10^(par("usr")[4] - 0.02 * diff(par("usr")[3:4]))
  text(mx, my, mstr, adj = c(0, 1), xpd = NA)
  ## Initial parameter estimates ...
  ## hacking to get alignment at "=" and at "e"
  pstr1 <- paste0(names(x$theta_mle), " = ")
  mat <- matrix(unlist(strsplit(sprintf("%0.3e", unlist(x$theta_mle)), "e")),
                nrow = 2)
  pstr2 <- paste0(mat[1, ], "e")
  pstr3 <- paste0(mat[2, ])
  px1 <- mx + max(strwidth(pstr1))
  px2 <- px1 + max(strwidth(pstr2))
  py <- 10^(par("usr")[3] + 0.02 * diff(par("usr")[3:4]))
  text(px1, py, paste(pstr1, collapse = "\n"), adj = c(1, 0), xpd = NA)
  text(px2, py, paste(pstr2, collapse = "\n"), adj = c(1, 0), xpd = NA)
  text(px2, py, paste(pstr3, collapse = "\n"), adj = c(0, 0), xpd = NA)
  par(op)
  invisible(NULL)
}

