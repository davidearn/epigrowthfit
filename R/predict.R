#' Compute predicted values
#'
#' Computes predicted values of interval incidence, cumulative incidence,
#' and the per capita growth rate, conditional on observed data and a fitted
#' nonlinear mixed effects model of epidemic growth.
#'
#' @param object
#'   An \code{"\link{egf}"} object.
#' @param what
#'   A \link{character} vector listing one or more variables for which
#'   predicted values are sought.
#' @param time
#'   A \link{numeric} vector supplying time points at which predicted
#'   values are sought. \link{Date} and \link{POSIXt} vectors are
#'   tolerated and coerced to numeric with \code{\link{julian}(time)}.
#'   When \link{time} is missing, the original time points stored in
#'   \code{object} are used.
#' @param window
#'   A \link{factor} of length \code{\link{length}(time)} grouping
#'   the elements of \code{time} by fitting window. Levels not
#'   belonging to \code{\link{levels}(object$frame$window)} are ignored.
#'   \code{window} is ignored altogether if \code{time} is missing.
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
#' In the result, \code{estimate[i]} can be interpreted as follows
#' assuming \code{log = FALSE}:
#' \describe{
#' \item{\code{interval}}{
#'   The expected number of cases observed from \code{time[i-1]}
#'   to \code{time[i]} in \code{window[i]} (interval incidence).
#' }
#' \item{\code{cumulative}}{
#'   The expected number of cases observed up to \code{time[i]}
#'   in \code{window[i]} (cumulative incidence).
#' }
#' \item{\code{rt}}{
#'   The predicted per capita growth rate at \code{time[i]}.
#'   This is computed exactly from the differential equation model
#'   associated with \code{object$model$curve}.
#' }
#' }
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
#'   Time series, from \code{\link{levels}(object$frame_windows$ts)}.
#' }
#' \item{window}{
#'   Fitting window, from \code{\link{levels}(object$frame_windows$window)}.
#' }
#' \item{time}{
#'   Time, after possible coercion to \link{numeric}.
#' }
#' \item{estimate}{
#'   Predicted value of \code{var} at \code{time} in \code{window},
#'   conditional on the mixed effects data \code{object$frame_parameters}
#'   and the fitted model \code{object$best}.
#' }
#' \item{se}{
#'   (If \code{log = TRUE} and \code{se = TRUE}.)
#'   Approximate (delta method) standard error on \code{estimate}.
#' }
#'
#' @export
#' @importFrom TMB MakeADFun sdreport
predict.egf <- function(object,
                        what = c("interval", "cumulative", "rt"),
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
  append <- eval_append(substitute(append), combined, baseenv(), .append = .append)
  do_day_of_week <- object$model$day_of_week > 0L

  start <- object$frame_windows$start
  end <- object$frame_windows$end
  day1 <- object$tmb_args$data$day1

  if (missing_time <- missing(time)) {
    len <- object$tmb_args$data$time_seg_len
    time <- object$tmb_args$data$time + rep.int(start, len)
    time_split <- split(time, rep.int(seq_along(len), len))
    subset <- rep_len(TRUE, length(len))
  } else {
    if (inherits(time, c("Date", "POSIXt"))) {
      time <- julian(time)
    } else {
      stopifnot(is.numeric(time))
    }
    stopifnot(
      length(time) > 0L,
      is.finite(time)
    )
    if (do_day_of_week) {
      stopifnot(all.equal(time, z <- round(time)))
      time <- z
    }
    stopifnot(
      is.factor(window),
      length(window) == length(time)
    )

    subset <- levels(object$frame$window) %in% levels(factor(window))
    window <- factor(window, levels = levels(object$frame$window)[subset])

    stop_if_not(
      nlevels(window) > 0L,
      m = "`window` must have at least one valid level."
    )
    len <- c(table(window))
    min_len <- 1L + as.integer("interval" %in% what)
    stop_if_not(
      len >= min_len,
      m = sprintf("`time` must have length %d or greater in each level of `window`.", min_len)
    )

    time_split <- split(time, window)
    t0 <- vapply(time_split, min, 0)
    t1 <- vapply(time_split, max, 0)
    start <- start[subset]
    end <- end[subset]
    day1 <- day1[subset]

    stop_if_not(
      t0 >= start,
      t1 <= end,
      m = "`time[i]` must not occur before (after) the start (end) of `window[i]`."
    )
    if (do_day_of_week) {
      check_ok_diff_time <- function(x) all(diff(x) == 1)
    } else {
      check_ok_diff_time <- function(x) all(diff(x) > 0)
    }
    stop_if_not(
      vapply(time_split, check_ok_diff_time, FALSE),
      m = wrap(
        "`time` must be increasing ",
        if (do_day_of_week) "with one day spacing " else "",
        "in each level of `window`."
      )
    )

    time <- unlist(time_split, FALSE, FALSE)
    if (do_day_of_week) {
      day1 <- as.integer((day1 + (t0 - start)) %% 7)
    }
  }

  object$tmb_args$data$flags$flag_predict <- 1L
  object$tmb_args$data$flag_what <- as.integer(eval(formals(predict.egf)$what) %in% what)
  object$tmb_args$data$subset <- which(subset) - 1L
  object$tmb_args$data$new_time <- time - rep.int(start, len)
  object$tmb_args$data$new_time_seg_len <- len
  object$tmb_args$data$new_day1 <- day1
  tmb_out <- do.call(MakeADFun, object$tmb_args)

  if (log && se) {
    tmb_out$fn(object$best[object$nonrandom])
    ssdr <- summary(sdreport(tmb_out), select = "report")
    index <- factor(rownames(ssdr), levels = sprintf("log_%s", what), labels = sprintf("log(%s)", what))
    report <- split(unname(ssdr[, "Estimate"]), index)
    report_se <- split(unname(ssdr[, "Std. Error"]), index)
  } else {
    report <- tmb_out$report(object$best)[sprintf("log_%s", what)]
    names(report) <- sprintf("log(%s)", what)
  }

  last <- cumsum(len)
  first <- c(0L, last[-length(last)]) + 1L
  x <- rep_len(NA_real_, length(time))
  ix <- list(
    interval = -first,
    cumulative = seq_along(time),
    rt = seq_along(time)
  )

  res <- data.frame(
    var = gl(length(report), length(time), labels = names(report)),
    ts = rep.int(object$frame_windows$ts[subset], len),
    window = rep.int(object$frame_windows$window[subset], len),
    time = time,
    estimate = unlist(Map(replace, list(x), ix[what], report), FALSE, FALSE)
  )
  if (log && se) {
    res$se <- unlist(Map(replace, list(x), ix[what], report_se), FALSE, FALSE)
  }
  if (!log) {
    res$estimate <- exp(res$estimate)
    levels(res$var) <- what
  }
  res <- data.frame(
    res,
    combined[rep.int(which(subset), len), append, drop = FALSE],
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  attr(res, "se") <- log && se
  class(res) <- c("egf_predict", "data.frame")
  res
}

#' Confidence intervals on predicted values
#'
#' Computes confidence intervals on predicted values of interval incidence,
#' cumulative incidence, and the per capita growth rate.
#'
#' @param object
#'   An \code{"\link[=predict.egf]{egf_predict}"} object.
#'   Must supply log scale predicted values and corresponding standard errors.
#' @param parm
#'   Unused argument included for generic consistency.
#' @param level
#'   A number in the interval (0,1) indicating a confidence level.
#' @param log
#'   A \link{logical} flag. If \code{FALSE}, then confidence intervals
#'   on inverse log-transformed predicted values are returned.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' Confidence limits on predicted values (log scale) are computed
#' as \code{estimate + c(-1, 1) * sqrt(q) * se},
#' with \code{estimate} and \code{se} as in \code{object} and
#' \code{q = \link{qchisq}(level, df = 1)}.
#'
#' @return
#' If \code{log = TRUE}, then \code{object} but with variable
#' \code{se} replaced with variables \code{lower} and \code{upper}
#' supplying confidence limits on log predicted values.
#'
#' Otherwise, the same result but with variables \code{estimate},
#' \code{lower}, and \code{upper} inverse log-transformed and
#' \code{\link{levels}(var)} modified accordingly.
#'
#' \code{level} is retained as an \link[=attributes]{attribute}.
#'
#' @export
confint.egf_predict <- function(object, parm, level = 0.95, log = TRUE, ...) {
  if (!isTRUE(attr(object, "se"))) {
    stop(wrap(
      "`object` must supply log scale predicted values ",
      "and corresponding standard errors. ",
      "Retry with `object = predict(., link = TRUE, log = TRUE)`."
    ))
  }
  stop_if_not_number_in_interval(level, 0, 1, "()")
  stop_if_not_true_false(log)

  s <- c("var", "ts", "window", "estimate", "se")
  res <- data.frame(
    object[s[1:4]],
    do_wald(estimate = object$estimate, se = object$se, level = level),
    object[-match(s, names(object), 0L)],
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  attr(res, "level") <- level
  if (log) {
    return(res)
  }
  elu <- c("estimate", "lower", "upper")
  res[elu] <- exp(res[elu])
  levels(res$var) <- sub("^log\\((.+)\\)$", "\\1", levels(res$var))
  res
}
