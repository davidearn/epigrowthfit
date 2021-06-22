#' Compute predicted values
#'
#' Computes predicted values of interval incidence, cumulative incidence,
#' and the per capita growth rate conditional on observed data and an estimated
#' model of epidemic growth.
#'
#' @param object
#'   An \code{"\link{egf}"} object.
#' @param what
#'   A \link{character} vector listing one or more variables for which
#'   predicted values are sought.
#' @param time
#'   A \link{numeric} or \link{Date} vector supplying time points at which
#'   predicted values are sought. Numeric \code{time} is assumed to measure
#'   time as a number of days since \emph{the end of} Date
#'   \code{origin} = \link{attr}(object$frame, "origin")}. Date \code{time}
#'   is coerced to numeric \code{\link{julian}(time, origin)}. When missing,
#'   \code{object$frame$time} is used.
#' @param window
#'   A \link{factor} of length \code{\link{length}(time)} such that
#'   \code{\link{split}(time, window)} splits \code{time} by fitting window.
#'   Levels not belonging to \code{\link{levels}(object$frame$window)}
#'   are ignored. When \code{time} is missing, \code{object$frame$window}
#'   is used.
#' @param append
#'   An expression indicating variables in the combined model frame
#'   (see \code{\link{make_combined}}) to be included with the result.
#'   The default (\code{\link{NULL}}) is to append nothing.
#' @param log
#'   A \link{logical} flag. If \code{FALSE},
#'   then inverse log-transformed predicted values are returned.
#' @param se
#'   A \link{logical} flag. If \code{log = TRUE} and \code{se = TRUE},
#'   then approximate (delta method) standard errors on predicted values
#'   are reported. Standard errors are required for subsequent use of
#'   \code{\link{confint.egf_predict}}.
#' @param .append
#'   A \link{character} vector listing variable names to be used
#'   (if non-\code{\link{NULL}}) in place of the result of evaluating
#'   \code{append}.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' In the result, \code{exp(estimate[i])} can be interpreted as follows,
#' assuming \code{log = TRUE}:
#' \describe{
#' \item{\code{log_int_inc}}{
#'   The number of cases predicted from \code{time[i-1]} to \code{time[i]}
#'   in \code{window[i]} (interval incidence).
#' }
#' \item{\code{log_cum_inc}}{
#'   The number of cases predicted from the start of \code{window[i]}
#'   to \code{time[i]} (cumulative incidence).
#' }
#' \item{\code{log_rt}{
#'   The per capita growth rate at \code{time[i]}. This is obtained
#'   exactly from the differential equation model associated with
#'   \code{object$model$curve}, unless \code{object$model$day_of_week > 0},
#'   in which case it is approximated by 7-point local linear regression
#'   of the log interval incidence time series.
#' }
#' }
#'
#' See topic \code{\link{nse}} for details on nonstandard evaluation
#' of \code{append}.
#'
#' @return
#' A \link[=data.frame]{data frame} inheriting from \link{class}
#' \code{"egf_predict"}, with variables:
#' \item{var}{
#'   Predicted variable, from \code{what}.
#' }
#' \item{ts}{
#'   Time series, from \code{\link{levels}(object$endpoints$ts)}.
#' }
#' \item{window}{
#'   Fitting window, from \code{\link{levels](object$endpoints$window)}.
#' }
#' \item{time}{
#'   Time, after possible coercion from \link{Date} to \link{numeric}.
#' }
#' \item{estimate}{
#'   Predicted value of \code{var} at \code{time} in \code{window},
#'   conditional on the mixed effects data found in \code{object$frame_par}
#'   and the fitted model specified by \code{object$best}.
#' }
#' \item{se}{
#'   (If \code{log = TRUE} and \code{se = TRUE}.)
#'   Approximate (delta method) standard error on \code{estimate}.
#' }
#' \code{\link{attr}(object$frame, "origin")} is retained
#' as an \link[=attributes]{attribute}.
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
  stop_if_not_true_false(log)
  stop_if_not_true_false(se)

  combined <- make_combined(object)
  append <- append_to_index(substitute(append), data = combined, enclos = parent.frame(),
                            .append = .append)
  endpoints <- object$endpoints
  origin <- attr(endpoints, "origin")
  day_of_week <- object$model$day_of_week
  min_window_len <- max(2L, 8L * (day_of_week > 0L) * ("log_rt" %in% what))

  if (no_time <- missing(time)) {
    time <- object$frame$time
    window <- object$frame$window
    subset <- seq_along(levels(window))
  } else {
    if (inherits(time, "Date")) {
      time <- julian(time, origin = origin)
    }
    stop_if_not(
      is.numeric(time),
      is.finite(time),
      m = "`time` must be a finite numeric or Date vector."
    )
    if (day_of_week > 0L) {
      stop_if_not(
        all.equal(time, z <- round(time)),
        m = "object$model$day_of_week > 0: `time` must be an integer or Date vector."
      )
      time <- z
    }
    stop_if_not(
      is.factor(window),
      length(window) == length(time),
      m = "`window` must be a factor of length `length(time)`."
    )
    wl <- levels(object$frame_ts$window)
    subset <- match(levels(factor(window)), wl, 0L)
    window <- factor(window, levels = wl[subset])
    stop_if_not(
      nlevels(window) > 0L,
      m = "`window` must have at least one valid level."
    )
  }
  stop_if_not(
    tabulate(window) > min_window_len,
    m = sprintf("`time` must have length %d or greater in each level of `window`.", min_window_len)
  )
  time_split <- split(time, window)
  starts <- endpoints$start[subset]
  ends <- endpoints$end[subset]
  if (no_time) {
    stop_if_not(
      vapply(time_split, min, 0) >= starts,
      vapply(time_split, max, 0) <= ends,
      m = "`time[i]` must not occur before (after) the start (end) of `window[i]`."
    )
    if (day_of_week > 0L) {
      stop_if_not(
        vapply(time_split, function(x) all(diff(x) == 1), FALSE),
        n = "day_of_week > 0: `time` must have 1-day spacing in each level of `window`."
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
        if (day_of_week > 0L) {
          return(seq.int(x0, x[length(x)]))
        }
        if (x[1L] > x0) {
          return(c(x0, x))
        }
        x
      }
      time_split <- Map(f, x = time_split, x0 = starts)
    }
    time <- unlist(time_split, FALSE, FALSE)
  }

  ## Additional data objects needed to run prediction code in C++ template
  lens <- lengths(time_split, use.names = FALSE)
  l <- list(
    what_flag = as.integer(c("log_int_inc", "log_cum_inc", "log_rt") %in% what),
    t_predict = time - rep.int(starts, lens),
    t_predict_seg_len = lens,
    day_of_week_on_day0_predict = object$tmb_args$data$weekday_on_day0[subset],
    Yo_predict = object$tmb_args$data$Yo[subset, , drop = FALSE],
    Xs_predict = object$tmb_args$data$Xs[subset, , drop = FALSE],
    Xd_predict = object$tmb_args$data$Xd[subset, , drop = FALSE],
    Z_predict = object$tmb_args$data$Z[subset, , drop = FALSE]
  )
  object$tmb_args$data <- c(object$tmb_args$data, l)
  object$tmb_args$data$predict_flag <- 1L
  tmb_out <- do.call(MakeADFun, object$tmb_args)

  if (se) {
    tmb_out$fn(object$best[object$nonrandom])
    r <- as.list(sdreport(tmb_out))[what]
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
