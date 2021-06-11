#' Compute predicted values
#'
#' Computes predicted values of interval incidence,
#' cumulative incidence, and the per capita growth rate
#' conditional on observed data and an estimated model.
#'
#' @param object
#'   An `"egf"` object returned by [egf()].
#' @param what
#'   A character vector listing one or more variables for which
#'   predicted values are sought.
#' @param time
#'   A numeric or Date vector supplying time points at which
#'   predicted values are sought. Numeric `time` is assumed
#'   to measure time as a number of days since _the end of_
#'   date `origin = attr(object$frame_ts, "origin")`.
#'   Date `time` is coerced to numeric `julian(time, origin)`.
#'   When missing, `object$frame_ts$time` is used.
#' @param window
#'   A factor of length `length(time)` such that
#'   `split(time, window)` splits `time` by fitting window.
#'   Levels not belonging to `levels(object$frame_ts$window)`
#'   are ignored.
#'   When `time` is missing, `object$frame_ts$window` is used.
#' @param log
#'   A logical scalar. If `FALSE`, then inverse log-transformed
#'   predicted values are returned.
#' @param se
#'   A logical scalar. If `log = TRUE` and `se = TRUE`, then
#'   approximate (delta method) standard errors on predicted values
#'   are reported. Note that standard errors are required for
#'   subsequent use of [confint.egf_predict()].
#' @inheritParams fitted.egf
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' In the returned data frame, assuming `log = FALSE`,
#' `estimate[i]` can be interpreted as follows:
#' \describe{
#' \item{`int_inc`}{
#'   The number of cases predicted from `time[i-1]` to `time[i]`
#'   in `window[i]` (interval incidence).
#' }
#' \item{`cum_inc`}{
#'   The number of cases predicted from the start of `window[i]`
#'   to `time[i]` (cumulative incidence).
#' }
#' \item{`rt`}{
#'   The per capita growth rate at `time[i]`. This is obtained
#'   exactly from the differential equation model associated with
#'   `object$curve`, unless `object$weekday > 0`, in which case
#'   it is approximated by 7-point local linear regression of the
#'   log interval incidence time series.
#' }
#' }
#'
#' See topic [`nse`] for details on nonstandard evaluation of `append`.
#'
#' @return
#' A data frame inheriting from class `"egf_predict"`, with variables:
#' \item{`var`}{
#'   Predicted variable, from `what`.
#' }
#' \item{`ts`}{
#'   Time series, from `levels(object$endpoints$ts)`.
#' }
#' \item{`window`}{
#'   Fitting window, from `levels(object$endpoints$window)`.
#' }
#' \item{`time`}{
#'   Time, after possible coercion from Date to numeric.
#' }
#' \item{`estimate`}{
#'   Predicted value of `var` at `time` in `window` (see Details),
#'   conditional on `object$frame_par` and `object$best`.
#' }
#' \item{`se`}{
#'   (If `log = TRUE` and `se = TRUE`.)
#'   Approximate (delta method) standard error on `estimate`.
#' }
#' `attr(object$endpoints, "origin")` is retained as an attribute.
#'
#' @export
#' @importFrom TMB MakeADFun sdreport
predict.egf <- function(object,
                        what = c("log_int_inc", "log_cum_inc", "log_rt"),
                        time,
                        window,
                        append = NULL,
                        log = TRUE,
                        se = FALSE,
                        .append = NULL,
                        ...) {
  what <- unique(match.arg(what, several.ok = TRUE))
  stop_if_not_true_false(se)

  combined <- make_combined(object)
  append <- append_to_index(substitute(append), combined, parent.frame(),
                            .append = .append)

  endpoints <- object$endpoints
  origin <- attr(endpoints, "origin")
  weekday <- object$weekday
  min_window_len <- max(2L, 8L * (weekday > 0L) * ("log_rt" %in% what))

  if (no_time <- missing(time)) {
    time <- object$frame_ts$time
    window <- object$frame_ts$window
    subset <- seq_along(levels(window))
  } else {
    if (inherits(time, "Date")) {
      time <- julian(time, origin = origin)
    }
    stop_if_not(
      is.numeric(time),
      m = "`time` must be a numeric or Date vector."
    )
    stop_if_not(
      !anyNA(time),
      m = "`time` must not have missing values."
    )
    if (weekday > 0L) {
      stop_if_not(
        all.equal(time, z <- round(time)),
        m = paste0(
          "weekday > 0: `time` must be an integer vector\n",
          "(in the sense of `all.equal(time, round(time))`)."
        )
      )
      time <- z
    }
    stop_if_not(
      is.factor(window),
      length(window) == length(time),
      m = "`window` must be a factor of length `length(time)`."
    )
    wl <- levels(object$frame_ts$window)
    subset <- match(levels(droplevels(window)), wl, 0L)
    window <- factor(window, levels = wl[subset])
    stop_if_not(
      nlevels(window) > 0L,
      m = "`window` must have at least one valid level."
    )
  }
  stop_if_not(
    tabulate(window) > min_window_len,
    m = sprintf("`time` must have length %d or greater\nin each level of `window`.", min_window_len)
  )
  time_split <- split(time, window)
  starts <- endpoints$start[subset]
  ends <- endpoints$end[subset]
  if (no_time) {
    stop_if_not(
      vapply(time_split, min, 0) >= starts,
      vapply(time_split, max, 0) <= ends,
      m = "`time[i]` must not occur before (after)\nthe start (end) of `window[i]`."
    )
    if (weekday > 0L) {
      stop_if_not(
        vapply(time_split, function(x) all(diff(x) == 1), FALSE),
        n = "weekday > 0: `time` must have 1-day spacing\nin each level of `window`."
      )
    } else {
      stop_if_not(
        vapply(time_split, function(x) all(diff(x) > 0), FALSE),
        m = "`time` must be increasing in each level of `window`."
      )
    }
    ## To compute cumulative incidence since the start of a window,
    ## the start of the window must be included as a time point
    if ("log_cum_inc" %in% what) {
      f <- function(x, x0) {
        if (weekday > 0L) {
          return(seq.int(x0, x[length(x)]))
        }
        if (x[1L] > x0) {
          return(c(x0, x))
        }
        x
      }
      time_split <- Map(f, x = time_split, x0 = starts)
    }
    time <- unlist(time_split, use.names = FALSE)
  }

  ## Additional data objects needed to run prediction code
  ## in C++ template
  lens <- lengths(time_split, use.names = FALSE)
  l <- list(
    what_flag = as.integer(c("log_int_inc", "log_cum_inc", "log_rt") %in% what),
    se_flag = as.integer(se),
    t_predict = time - rep.int(starts, lens),
    t_predict_seg_len = lens,
    weekday_on_day0_predict = object$tmb_args$data$weekday_on_day0[subset],
    Xd_predict = object$tmb_args$data$Xd[subset, , drop = FALSE],
    Xs_predict = object$tmb_args$data$Xs[subset, , drop = FALSE],
    Z_predict = object$tmb_args$data$Z[subset, , drop = FALSE],
    Yo_predict = object$tmb_args$data$Yo[subset, , drop = FALSE]
  )
  object$tmb_args$data$predict_flag <- 1L
  object$tmb_args$data <- c(object$tmb_args$data, l)
  tmb_out <- do.call(MakeADFun, c(object$tmb_args))

  if (se) {
    tmb_out$fn(object$best[object$nonrandom])
    r <- split_sdreport(sdreport(tmb_out))[what]
  } else {
    r <- tmb_out$report(object$best)[what]
    r <- lapply(r, function(x) list(estimate = x))
  }

  out <- list()
  x <- rep_len(NA_real_, length(time))
  el <- data.frame(
    var = rep_len(factor("foo"), length(time)),
    ts = rep.int(endpoints$ts[subset], lens),
    window = rep.int(endpoints$window[subset], lens),
    time = time,
    estimate = x,
    se = x,
    combined[rep.int(subset, lens), append, drop = FALSE]
  )
  lasts <- cumsum(lens)
  firsts <- c(0L, lasts[-length(lasts)]) + 1L
  edges <- unlist(Map(function(a, b) c(a + 0:3, b - 2:0), a = firsts, b = lasts))

  if ("log_int_inc" %in% what) {
    levels(el$var) <- "log_int_inc"
    el$estimate <- replace(x, -firsts, r$log_int_inc$estimate)
    if (se) {
      el$se <- replace(x, -firsts, r$log_int_inc$se)
    }
    out$log_int_inc <- el
  }
  if ("log_cum_inc" %in% what) {
    levels(el$var) <- "log_cum_inc"
    el$estimate <- replace(x, -firsts, r$log_cum_inc$estimate)
    if (se) {
      el$se <- replace(c, -firsts, r$log_cum_inc$se)
    }
    out$log_cum_inc <- el
  }
  if ("log_rt" %in% what) {
    levels(el$var) <- "log_rt"
    if (weekday > 0L) {
      el$estimate <- replace(x, -edges, r$log_rt$estimate)
      if (se) {
        el$se <- replace(x, -edges, r$log_rt$se)
      }
    } else {
      el$estimate <- r$log_rt$estimate
      if (se) {
        el$se <- r$log_rt$se
      }
    }
    out$log_rt <- el
  }
  out <- do.call(rbind, out)
  row.names(out) <- NULL
  attr(out, "origin") <- origin
  attr(out, "se") <- se
  class(out) <- c("egf_predict", "data.frame")
  out
}

#' Confidence intervals on predicted values
#'
#' Computes confidence intervals on predicted values of
#' log cumulative incidence, log interval incidence, and
#' log per capita growth rate.
#'
#' @param object
#'   An `"egf_predict"` object returned by [predict.egf()].
#'   Must supply standard errors on predicted values.
#' @param log
#'   A logical scalar. If `FALSE`, then confidence intervals
#'   on inverse log-transformed predicted values are returned.
#' @inheritParams confint.egf
#'
#' @details
#' Confidence limits on predicted values (log scale) are computed
#' as `estimate + c(-1, 1) * sqrt(q) * se`, with `estimate` and
#' `se` obtained from `object` and `q = qchisq(level, df = 1)`.
#'
#' @return
#' If `log = TRUE`, then `object` but with variable `se` replaced
#' with variables `lower` and `upper` supplying confidence limits
#' on log predicted values.
#'
#' Otherwise, the same object but with variables `estimate`, `lower`,
#' and `upper` inverse log transformed and prefix `"log_"` stripped
#' from `levels(var)`.
#'
#' `level` is retained as an attribute.
#'
#' @export
confint.egf_predict <- function(object, parm, level = 0.95, log = TRUE, ...) {
  stop_if_not(
    attr(object, "se"),
    m = paste0(
      "`object` must supply standard errors on predicted values.\n",
      "Repeat `predict()` with `log = TRUE` and `se = TRUE`."
    )
  )
  stop_if_not_number_in_interval(level, 0, 1, "()")
  stop_if_not_true_false(log)

  d <- data.frame(
    object[c("var", "ts", "window", "estimate")],
    do_wald(estimate = object$estimate, se = object$se, level = level),
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  attr(d, "level") <- level
  if (log) {
    return(d)
  }
  s_elu <- c("estimate", "lower", "upper")
  d[s_elu] <- exp(d[s_elu])
  levels(d$par) <- sub("^log_", "", (levels(d$par)))
  d
}
