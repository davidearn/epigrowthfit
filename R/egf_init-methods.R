#' Methods for class "egf_init"
#'
#' @description
#' Methods for printing and plotting "egf_init" objects
#' returned by [egf_init()].
#'
#' @param x An "egf_init" object.
#' @param inc One of `"interval"` and `"cumulative"`,
#'   indicating whether to plot `x$cases` (interval incidence)
#'   or `cumsum(x$cases)` (cumulative incidence).
#' @param tol A non-negative number used only if
#'   `inc = "interval"`. `x$cases[i]` is plotted
#'   in light blue if `diff(x$time)[i] < (1-tol)*m`,
#'   in dark blue if `diff(x$time)[i] > (1+tol)*m`,
#'   and in grey otherwise, where `m = median(diff(x$time))`.
#'   Assign `Inf` to ensure that everything is grey.
#' @param ... Unused optional arguments.
#'
#' @details
#' If `x$time` is not equally spaced, then `inc = "interval"`
#' should be used with caution. `x$cases[i]` is the number
#' of cases observed between `x$time[i]` and `x$time[i+1]`,
#' hence `x$cases` roughly scales with `diff(x$time)`.
#' Argument `tol` can be used with `inc = "interval"`
#' in order to highlight outliers in `diff(x$time)`.
#'
#' @name egf_init-methods
NULL

#' @rdname egf_init-methods
#' @export
print.egf_init <- function(x, ...) {
  if (!inherits(x, "egf_init")) {
    stop("`x` must be an \"egf_init\" object.")
  }

  cstr <- switch(x$curve,
    exponential = "an exponential model",
    logistic    = "a logistic model",
    richards    = "a Richards model"
  )
  bstr <- if (x$include_baseline) "with a linear baseline and" else "with"
  dstr <- switch(x$distr,
    poisson = "Poisson-distributed observations.",
    nbinom  = "negative binomial observations."
  )
  f <- x$first
  l <- x$last
  cat("Pass this \"egf_init\" object to `egf()` to fit", cstr, "\n")
  cat(bstr, dstr, "\n")
  cat("\n")
  cat("Fitting window:\n")
  cat("\n")
  cat("index   ", f, ":", l, "\n", sep = "")
  cat(" time   [", x$time[f+1], ", ", x$time[l+1], "] days since ", as.character(x$date[1]), "\n", sep = "")
  cat("cases   ", sum(x$cases[f:l]), " of ", sum(x$cases), "\n", sep = "")
  cat("\n")
  cat("Initial parameter estimates:\n")
  cat("\n")
  print(unlist(x$theta0))
  invisible(x)
}

#' @rdname egf_init-methods
#' @export
#' @import graphics
#' @importFrom stats median
plot.egf_init <- function(x, inc = "interval", tol = 0.025, ...) {
  if (!inherits(x, "egf_init")) {
    stop("`x` must be an \"egf_init\" object.")
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

  op <- par(mar = c(5, 4, 4, 8) + 0.1, las = 1, mgp = c(3, 0.7, 0))
  data <- data.frame(time = x$time[-1], cases = x$cases)
  xlab <- paste("days since", as.character(x$date[1]))
  if (inc == "interval") {
    dtime <- diff(x$time)
    m <- median(dtime)
    dtime_min <- (1 - tol) * m
    dtime_max <- (1 + tol) * m
    bg_enum <- 1 + 1 * (dtime < dtime_min) + 2 * (dtime > dtime_max)
    bg <- c("#DDDDDD", "#66CCEE", "#4477AA")[bg_enum]
    plot(cases + 0.1 ~ time, data = data, xlab = xlab,
         log = "y", pch = 21, bg = bg)
  } else if (inc == "cumulative") {
    plot(cumsum(cases) + 0.1 ~ time, data = data, xlab = xlab,
         log = "y", pch = 21, bg = "#DDDDDD")
  }
  f <- x$first
  l <- x$last
  abline(v = x$time[c(f,l)+1], lty = 2, col = "#555555")
  axis(side = 3, at = x$time[c(f,l)+1], labels = c(f,l),
       tick = FALSE, mgp = c(3, 0.1, 0))
  mtext("index", side = 3, line = 2)
  ## Model
  mstr <- paste0(x$curve,
                 if (x$include_baseline) "\nbaseline" else "",
                 "\n", x$distr)
  mx <- par("usr")[2] + 0.02 * diff(par("usr")[1:2])
  my <- 10^(par("usr")[4] - 0.02 * diff(par("usr")[3:4]))
  text(mx, my, mstr, adj = c(0, 1), xpd = NA)
  ## Initial parameter estimates ...
  ## hacking to get alignment at "=" and at "e"
  pstr1 <- paste0(names(x$theta0), " = ")
  mat <- matrix(unlist(strsplit(sprintf("%0.3e", unlist(x$theta0)), "e")),
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

