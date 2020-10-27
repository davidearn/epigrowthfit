#' Methods for class "egf_init"
#'
#' @description
#' Methods for "egf_init" objects returned by [egf_init()].
#'
#' @param x,object An "egf_init" object.
#' @param log A logical scalar. For the `coef` method, if `TRUE`,
#'   then parameter estimates are log-transformed. For the plot
#'   method, if `TRUE`, then incidence is displayed on a
#'   logarithmic scale.
#' @param time A numeric vector listing increasing time points,
#'   expressed as numbers of days since `object$date[1]`.
#'   Missing values are not tolerated.
#' @param inc One of `"cumulative"` and `"interval"`,
#'   indicating a type of incidence to plot.
#' @param xty One of `"Date"` and `"numeric"`, indicating
#'   a type of horizontal axis.
#' @param add A logical scalar. If `TRUE`, then the fitting window
#'   and predicted curve are added to the current graphics device,
#'   presumably started by a previous invocation of the plot method.
#' @param annotate A logical scalar. If `TRUE`, then a legend and
#'   a list of parameter estimates are displayed in the right margin.
#' @param tol A non-negative number, used only if `inc = "interval"`.
#'   `cases[i]` is highlighted according to `point_style_short` if
#'   `diff(time)[i] < (1-tol)*m` and according to `point_style_long`
#'   if `diff(time)[i] > (1+tol)*m`, where `m = median(diff(time))`.
#'   In both cases, the value of `diff(time)[i]` is printed above
#'   the point. Assign 0 to ensure that all deviations from `m` are
#'   highlighted. Assign `Inf` to disable highlighting.
#' @param date_style A named list of arguments to [daxis()]
#'   (a subset of `tcl`, `line`, `col.axis`, and `cex.axis`),
#'   affecting the appearance of the bottom axis if `xty = "Date"`.
#' @param window_style A named list of arguments to
#'   [`polygon()`][graphics::polygon()], affecting the appearance
#'   of the fitting window. Currently, only `col` and `border` are
#'   implemented.
#' @param line_style A named list of arguments to
#'   [`lines()`][graphics::lines()], affecting the appearance
#'   of the predicted curve. Currently, only `lty`, `lwd`, and `col`
#'   are implemented.
#' @param point_style_main A named list of arguments to
#'   [`points()`][graphics::points()], affecting the appearance
#'   of the observed data. Currently, only `pch`, `col`, `bg` and `cex`
#'   are implemented.
#' @param point_style_short,point_style_long
#'   Alternatives to `point_style_main` used to highlight
#'   certain points when `inc = "interval"`. See argument `tol`.
#' @param text_style A named list of arguments to
#'   [`text()`][graphics::text()], affecting the appearance
#'   of text printed above highlighted points. See argument `tol`.
#'   Currently, only `pos`, `offset`, `col`, `cex`, and `font`
#'   are implemented.
#' @param ... Optional arguments. Used only by the `plot` method
#'   to specify graphical parameters. Currently, only `xlim`,
#'   `ylim`, `xlab`, `ylab`, and `main` are implemented. Further
#'   arguments are ignored.
#'   See [`plot()`][graphics::plot()] and [`par()`][graphics::par()]
#'   for a catalogue of graphical parameters.
#'
#' @return
#' The `print` method returns `x` invisibly.
#'
#' The `coef` method returns `object$theta_init` if `log = FALSE`
#' and `object$log_theta_init` if `log = TRUE`.
#'
#' The `predict` method returns a list with numeric elements:
#'
#' \describe{
#'   \item{`time`}{Matches argument.}
#'   \item{`refdate`}{Matches `object$date[1]`.}
#'   \item{`cum_inc`}{Expected cumulative incidence at time points
#'     `time`, conditional on parameter vector `object$theta_init`.
#'     Equal to `object$eval_cum_inc(time)`.
#'   }
#'   \item{`int_inc`}{Expected interval incidence given interval
#'     endpoints `time`, conditional on parameter vector
#'     `object$theta_init`. Equal to `diff(object$eval_cum_inc(time))`
#'     if `length(time) >= 2` and omitted otherwise.
#'   }
#' }
#'
#' The `plot` method returns `NULL` invisibly.
#'
#' @details
#' ## Plot elements
#'
#' If `xty = "Date"`, then the bottom axis displays the dates specified
#' by `x$date`. If `xty = "numeric"`, then the bottom axis displays the
#' number of days since `x$date[1]`. Regardless of `xty`, numeric
#' coordinates are used, hence the left and right boundaries of the plot
#' region are specified by `range(x$time)`.
#'
#' The left axis measures interval or cumulative incidence (depending
#' on `inc`). When incidence is displayed on a logarithmic scale, zeros
#' are plotted as a positive number less than 1. They are therefore
#' distinguished from nonzero counts, which are always at least 1.
#'
#' Observed data, specified by `x$date` and either `x$cases` or
#' `cumsum(x$cases)` (depending on `inc`), are plotted as points.
#' `cases[i]` gives the number of cases observed between `date[i]`
#' and `date[i+1]`, and `cumsum(cases)[i]` the number observed
#' between `date[1]` and `date[i+1]`. Both are plotted at `time[i+1]`.
#'
#' The fitting window, specified by indices `x$first` and `x$last`,
#' is displayed as a shaded rectangle behind the other plot elements.
#' The left and right boundaries occur at `time[first]` and
#' `time[last+1]`. (`cases[first]` is a count from `time[first]`
#' to `time[first+1]`, and `cases[last]` is a count from `time[last]`
#' to `time[last+1]`. Hence the fitting window starts at `time[first]`
#' and ends at `time[last+1]`.)
#'
#' The incidence curve predicted by initial parameter estimates
#' `x$theta_init` is displayed as a line supported on grid points
#' `wgrid = seq(time[first], time[last+1], by)`, where `by = 1`
#' for cumulative incidence and `by = median(diff(time))` for
#' interval incidence (ensuring that the interval incidence
#' curve has the correct scale; see below). The predicted curve
#' is obtained with `predict(x, wgrid)`.
#'
#' Careful interpretation of interval incidence is required if
#' the plotted time series is not equally spaced, because `cases`
#' roughly scales with `diff(time)`. That is, certain observations
#' may vary from the predicted curve not due to chance, but because
#' they represent a count over fewer or more days than the typical
#' observation interval, namely `median(diff(time))`. Observations
#' for which `diff(time)` differs from `median(diff(time))` are
#' highlighted according to argument `tol` and labeled with the
#' value of `diff(time)`.
#'
#' If `annotate = TRUE` and `add = FALSE`, then a legend and the
#' initial parameter estimates `x$theta_init` are displayed in the
#' right margin.
#'
#' @seealso [egf_init()]
#'
#' @name egf_init-methods
NULL

#' @rdname egf_init-methods
#' @export
print.egf_init <- function(x, ...) {
  cstr <- switch(x$curve,
    exponential = "an exponential model",
    logistic    = "a logistic model",
    richards    = "a Richards model"
  )
  bstr <- if (x$include_baseline) {
    "with a linear baseline and"
  } else {
    "with"
  }
  dstr <- switch(x$distr,
    poisson = "Poisson-distributed observations.",
    nbinom  = "negative binomial observations."
  )
  uvec <- c(r = "per day", thalf = "days", b = "per day")
  uvec <- uvec[names(uvec) %in% names(x$theta_init)]
  cat("Pass this \"egf_init\" object to `egf()` to fit", cstr, "\n")
  cat(bstr, dstr, "\n")
  cat("\n")
  cat("Fitting window:\n")
  cat("\n")
  cat("index   ", x$first, ":", x$last, "\n", sep = "")
  cat(" date   (", as.character(x$date[x$first]), ", ", as.character(x$date[x$last+1]), "]\n", sep = "")
  cat("cases   ", sum(x$cases[x$first:x$last]), " of ", sum(x$cases), "\n", sep = "")
  cat("\n")
  cat("Initial parameter estimates:\n")
  cat("\n")
  print(x$theta_init)
  cat("\n")
  cat("Units:\n")
  cat("\n")
  print(uvec, quote = FALSE)
  invisible(x)
}

#' @rdname egf_init-methods
#' @export
coef.egf_init <- function(object, log = FALSE, ...) {
  if (!is.logical(log) || length(log) != 1 || is.na(log)) {
    stop("`log` must be `TRUE` or `FALSE`.")
  }

  if (log) {
    object$log_theta_init
  } else {
    object$theta_init
  }
}

#' @rdname egf_init-methods
#' @export
predict.egf_init <- function(object, time = object$time, ...) {
  if (!is.numeric(time) || length(time) == 0) {
    stop("`time` must be numeric and have nonzero length.")
  } else if (anyNA(time)) {
    stop("`time` must not have missing values.")
  } else if (!all(diff(time) > 0)) {
    stop("`time` must be increasing.")
  } else if (any(time < object$time[object$first] | time > object$time[object$last+1])) {
    warning("There are elements of `time` outside of the fitting window.",
             call. = FALSE)
  }

  out <- list(
    time = time,
    refdate = object$date[1],
    cum_inc = object$eval_cum_inc(time)
  )
  if (length(time) > 1) {
    out$int_inc = diff(out$cum_inc)
  }
  out
}

#' @rdname egf_init-methods
#' @export
#' @import graphics
#' @importFrom stats median
plot.egf_init <- function(x, inc = "interval", xty = "Date", log = TRUE,
                          add = FALSE, annotate = TRUE, tol = 0,
                          date_style = list(tcl = -0.2, line = c(0.05, 1), col.axis = c("black", "black"), cex.axis = c(0.7, 0.85)),
                          window_style = list(col = "#DDCC7740", border = NA),
                          line_style = list(lty = 1, lwd = 3, col = "#44AA99"),
                          point_style_main = list(pch = 21, col = "#BBBBBB", bg = "#DDDDDD", cex = 1),
                          point_style_short = list(pch = 1, col = "#882255", bg = NA, cex = 1),
                          point_style_long = list(pch = 16, col = "#882255", bg = NA, cex = 1),
                          text_style = list(pos = 3, offset = 0.3, col = "#BBBBBB", cex = 0.7, font = 2),
                          ...) {
  if (!is.character(inc) || length(inc) != 1 ||
      !inc %in% c("interval", "cumulative")) {
    stop("`inc` must be one of \"interval\", \"cumulative\".")
  }
  if (!is.character(xty) || length(xty) != 1 ||
      !tolower(xty) %in% c("date", "numeric")) {
    stop("`xty` must be one of \"Date\", \"numeric\".")
  }
  if (!is.logical(log) || length(log) != 1 || is.na(log)) {
    stop("`log` must be `TRUE` or `FALSE`.")
  }
  if (!is.logical(add) || length(add) != 1 || is.na(add)) {
    stop("`add` must be `TRUE` or `FALSE`.")
  }
  if (!is.logical(annotate) || length(annotate) != 1 || is.na(annotate)) {
    stop("`annotate` must be `TRUE` or `FALSE`.")
  }
  if (inc == "interval") {
    if (!is.numeric(tol) || length(tol) != 1 || !isTRUE(tol >= 0)) {
      stop("`tol` must be a non-negative number.")
    }
  }
  ss <- grep("_style", names(formals(plot.egf_init)), value = TRUE)
  for (s in ss) {
    l <- get(s)
    if (!is.list(l)) {
      stop("All \"_style\" arguments must be lists.")
    }
  }


  ### SET UP ###########################################################

  ## Optional graphical parameters
  dots <- list(...)

  ## Observed data
  data <- data.frame(
    time = x$time[-1],
    cum_inc = cumsum(x$cases),
    int_inc = x$cases,
    dt = diff(x$time)
  )
  m <- median(data$dt)

  ## Predicted curve
  wleft <- x$time[x$first]
  wright <- x$time[x$last+1]
  wgrid <- seq(wleft, wright, by = if (inc == "interval") m else 1)
  wpred <- predict(x, wgrid)[c("time", "cum_inc", "int_inc")]
  wpred$int_inc <- c(NA, wpred$int_inc)

  ## A way to mostly avoid conditional `if (inc = ...) ... else ...`
  varname <- substr(inc, start = 1, stop = 3) # first three characters
  varname <- paste0(varname, "_inc")
  formula <- as.formula(paste(varname, "~ time"))

  ## A way to include zeros on a logarithmic scale
  ymax <- max(data[[varname]])
  zero <- if (log) ymax^-0.04 else 0
  if (log) {
    data[[varname]][data[[varname]] == 0] <- zero
  }

  ## Titles
  if ("xlab" %in% names(dots)) {
    xlab <- dots$xlab
  } else if (tolower(xty) == "date") {
    xlab <- "date"
  } else if (xty == "numeric") {
    xlab <- paste("days since", as.character(x$date[1]))
  }
  if ("ylab" %in% names(dots)) {
    ylab <- dots$ylab
  } else {
    ylab <- paste(inc, "incidence")
  }
  cstr <- x$curve
  substr(cstr, 1, 1) <- toupper(substr(cstr, 1, 1)) # capitalize first letter
  if ("main" %in% names(dots)) {
    main <- dots$main
  } else {
    main <- paste(cstr, "model of", inc, "incidence\n(initial guess)")
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
    ymin <- zero
    ymax <- if (log) ymax^1.04 else ymax * 1.04
    ylim <- c(ymin, ymax)
  }

  ## Styles
  for (s in ss) {
    l1 <- eval(formals(plot.egf_init)[[s]]) # default style
    l2 <- get(s) # passed style
    inter <- intersect(names(l1), names(l2))
    l1[inter] <- l2[inter]
    assign(s, l1)
  }

  ## Style for each point
  ## (for interval incidence, style depends on observation interval)
  dt_min <- (1 - tol) * m
  dt_max <- (1 + tol) * m
  dt_enum <- 1 + (inc == "interval") * (1 * (data$dt < dt_min) + 2 * (data$dt > dt_max))
  pss <- c("main", "short", "long")
  data$point_style <- pss[dt_enum]


  ### PLOT #############################################################

  if (add) {
    sp <- get("par", envir = .egf_env)
    sp$yaxp <- NULL
    op <- par(sp)
  } else {
    op <- par(mar = c(4, 5, 2.7, 0.5 + 6 * annotate) + 0.1)
    plot.new()
    plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i",
                log = if (log) "y" else "")
  }

  ## Fitting window
  l <- list(
    x = c(wleft, wright, wright, wleft),
    y = ylim[c(1, 1, 2, 2)]
  )
  do.call(polygon, c(l, window_style))
  windex <- (data$time >= wleft - 4 & data$time <= wright + 4)

  if (!add) {
    ## Box
    box(bty = "l")

    ## Axis (x)
    if (tolower(xty) == "date") {
      l <- list(
        left = par("usr")[1],
        right = par("usr")[2],
        refdate = x$date[1]
      )
      do.call(daxis, c(l, date_style))
    } else if (xty == "numeric") {
      axis(side = 1, mgp = c(3, 0.7, 0), cex.axis = 0.85)
    }

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
    axis(side = 2, at = yax_at, labels = yax_labels, las = 1,
         mgp = c(3, 0.7, 0), cex.axis = if (digits > 1) 0.65 else 0.85)
  }

  ## Observed data
  for (ps in pss) {
    l <- list(
      formula = formula,
      data = data,
      subset = (data$point_style == ps) & (!add | windex),
      xpd = !any(c("xlim", "ylim") %in% names(dots))
    )
    do.call(points, c(l, get(paste0("point_style_", ps))))
  }

  ## Annotation above exceptional points
  if (any(dt_enum != 1)) {
    l <- list(
      formula = int_inc ~ time,
      data = data,
      labels = data$dt,
      subset = (dt_enum != 1) & (if (add) windex else TRUE)
    )
    do.call(text, c(l, text_style))
  }

  ## Predicted curve
  l <- list(
    formula = formula,
    data = wpred
  )
  do.call(lines, c(l, line_style))

  if (!add) {
    ## Titles
    title(xlab = xlab, line = 3)
    title(ylab = ylab, line = 4)
    title(main = main, line = 1, cex.main = 0.9)

    if (annotate) {
      ## Parameter estimates
      pstr1 <- paste0(names(x$theta_init), " = ")
      pstr2 <- round(x$theta_init, digits = 4)
      if ("K" %in% names(x$theta_init)) {
        pstr2[["K"]] <- round(pstr2[["K"]])
      }
      px <- par("usr")[2] + 0.02 * diff(par("usr")[1:2])
      px <- px + max(strwidth(pstr1, cex = 0.7))
      py <- par("usr")[3] + 0.02 * diff(par("usr")[3:4])
      if (log) {
        py <- 10^py
      }
      text(px, py, paste(pstr1, collapse = "\n"),
           adj = c(1, 0), xpd = NA, cex = 0.7)
      text(px, py, paste(pstr2, collapse = "\n"),
           adj = c(0, 0), xpd = NA, cex = 0.7)

      ## Legend (beware: some ugly hacks here)
      lx <- par("usr")[2] + 0.02 * diff(par("usr")[1:2])
      ly <- par("usr")[4] - 0.02 * diff(par("usr")[3:4])
      if (log) {
        ly <- 10^ly
      }
      if (inc == "cumulative") {
        lstr <- c("obs", NA, NA, "pred")
        index <- c(TRUE, FALSE, FALSE, TRUE)
      } else if (inc == "interval") {
        lstr1 <- paste0("'", c("obs,", "obs,", "obs,", "pred,"), "'")
        cond <- (all(data$dt[dt_enum == 1] == m))
        rel <- c(if(cond) "=" else "~", "<", ">", "=")
        mstr <- paste0(m, " day", if (m > 1) "s" else "")
        lstr2 <- paste0("'t ", rel, " ", mstr, "'")
        lstr <- parse(text = paste(lstr1, "~ Delta *", lstr2))
        index <- c(TRUE, any(dt_enum == 2), any(dt_enum == 3), TRUE)
      }
      legend(x = lx, y = ly,
             xpd = NA, bty = "n", cex = 0.7, seg.len = 1,
             legend = lstr[index],
             pch = c(point_style_main$pch, point_style_short$pch, point_style_long$pch, NA)[index],
             pt.bg = c(point_style_main$bg, point_style_short$bg, point_style_long$bg, NA)[index],
             lty = c(NA, NA, NA, line_style$lty)[index],
             lwd = c(NA, NA, NA, line_style$lwd)[index],
             col = c(point_style_main$col, point_style_short$col, point_style_long$col, line_style$col)[index])
    }
  }

  if (!add) {
    assign("par", par(no.readonly = TRUE), envir = .egf_env)
  }
  par(op)
  invisible(NULL)
}
