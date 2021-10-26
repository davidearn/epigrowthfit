#' Compute predicted values
#'
#' Computes predicted values of interval incidence, cumulative incidence,
#' and the per capita growth rate, conditional on observed data and a fitted
#' nonlinear mixed effects model of epidemic growth.
#'
#' @param object
#'   An \code{"\link{egf}"} object.
#' @param what
#'   A character vector listing one or more variables for which
#'   predicted values are sought.
#' @param time
#'   A numeric vector supplying time points at which predicted
#'   values are sought. \link{Date} and \link{POSIXt} vectors are
#'   tolerated and coerced to numeric with \code{\link{julian}(time)}.
#'   When \link{time} is missing, time points stored in
#'   \code{\link[=model.frame.egf]{model.frame}(object)} are reused,
#'   and \code{window} is ignored.
#' @param window
#'   A factor of length \code{length(time)} grouping the elements
#'   of \code{time} by fitting window.
#'   Levels not found in
#'   \code{levels(\link[=model.frame.egf]{model.frame}(object)$window)}
#'   are ignored.
#' @param log
#'   A logical flag. If \code{FALSE},
#'   then inverse log-transformed predicted values are returned.
#' @param se
#'   A logical flag. If \code{se = TRUE} and \code{log = TRUE},
#'   then approximate delta method standard errors on predicted values
#'   are reported.
#'   Standard errors are required for subsequent use
#'   of \code{\link{confint.egf_predict}}.
#'   Setting \code{se = TRUE} and \code{log = FALSE} is an error,
#'   as standard errors are not available for inverse log-transformed
#'   predicted values.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' In the result, \code{estimate[i]} can be interpreted as follows,
#' for \code{log = FALSE}:
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
#'
#' @return
#' A data frame inheriting from class \code{"egf_predict"}, with variables:
#' \item{var}{
#'   Predicted variable, from \code{what}.
#' }
#' \item{ts}{
#'   Time series, from
#'   \code{levels(\link[=model.frame.egf]{model.frame}(object)$ts)}.
#' }
#' \item{window}{
#'   Fitting window, from
#'   \code{levels(\link[=model.frame.egf]{model.frame}(object)$window)}.
#' }
#' \item{time}{
#'   Time, after possible coercion to numeric.
#' }
#' \item{estimate}{
#'   Predicted value of \code{var} at \code{time} in \code{window},
#'   conditional on the mixed effects data and fitted model.
#' }
#' \item{se}{
#'   Approximate delta method standard error on \code{estimate}.
#'   Absent except for calls matching \code{predict(log = TRUE, se = TRUE)}.
#' }
#'
#' @examples
#' object <- egf_cache("egf-1.rds")
#' zz <- egf_cache("predict-egf-1.rds", predict(object, se = TRUE))
#' str(zz)
#'
#' @export
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats model.frame
predict.egf <- function(object,
                        what = c("interval", "cumulative", "rt"),
                        time,
                        window,
                        log = TRUE,
                        se = FALSE,
                        ...) {
  what <- unique(match.arg(what, several.ok = TRUE))
  stopifnot(
    is_true_or_false(log),
    is_true_or_false(se)
  )
  if (se && !log) {
    stop(wrap(
      "Standard errors are not available for inverse log-transformed ",
      "predicted values."
    ))
  }

  frame_windows <- model.frame(object, "windows")
  start <- frame_windows$start
  end   <- frame_windows$end
  day1  <- object$tmb_out$env$data$day1
  do_day_of_week <- object$model$day_of_week > 0L

  if (missing_time <- missing(time)) {
    len  <- object$tmb_out$env$data$time_seg_len
    subset <- seq_along(len)
    time <- object$tmb_out$env$data$time + rep.int(start, len)
    time_split <- split(time, rep.int(subset, len))
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

    subset <- which(levels(frame_windows$window) %in% levels(factor(window)))
    window <- factor(window, levels = levels(frame_windows$window)[subset])

    if (nlevels(window) == 0L) {
      stop("'window' must have at least one valid level.")
    }
    len <- c(table(window))
    min_len <- 1L + as.integer("interval" %in% what)
    if (any(len < min_len)) {
      stop(sprintf("'time' must have length %d or greater in each level of 'window'.", min_len))
    }

    time_split <- split(time, window)
    t0 <- vapply(time_split, min, 0)
    t1 <- vapply(time_split, max, 0)

    start <- start[subset]
    end <- end[subset]
    day1 <- day1[subset]

    if (any(t0 < start | t1 > end)) {
      stop("'time[i]' must not occur before (after) the start (end) of 'window[i]'.")
    }
    if (do_day_of_week) {
      check_ok_diff_time <- function(x) all(diff(x) == 1)
    } else {
      check_ok_diff_time <- function(x) all(diff(x) > 0)
    }
    if (!all(vapply(time_split, check_ok_diff_time, FALSE))) {
      stop(wrap(
        "'time' must be increasing ",
        if (do_day_of_week) "with one day spacing " else "",
        "in each level of 'window'."
      ))
    }

    time <- unlist1(time_split)
    if (do_day_of_week) {
      day1 <- as.integer((day1 + (t0 - start)) %% 7)
    }
  }

  tmb_args <- egf_tmb_remake_args(object$tmb_out, par = object$best)
  tmb_args$data$flags$predict <- 1L
  tmb_args$data$what <- as.integer(eval(formals(predict.egf)$what) %in% what)
  tmb_args$data$subset <- subset - 1L
  tmb_args$data$new_time <- time - rep.int(start, len)
  tmb_args$data$new_time_seg_len <- len
  tmb_args$data$new_day1 <- day1
  tmb_out_retape <- do.call(MakeADFun, tmb_args)

  if (se) {
    sdr <- sdreport(tmb_out_retape, par.fixed = object$best[!object$random], getReportCovariance = FALSE)
    ssdr <- summary(sdr, select = "report")
    index <- factor(rownames(ssdr), levels = sprintf("log_%s", what), labels = sprintf("log(%s)", what))
    report <- split(unname(ssdr[, "Estimate"]), index)
    report_se <- split(unname(ssdr[, "Std. Error"]), index)
  } else {
    report <- tmb_out_retape$report(object$best)[sprintf("log_%s", what)]
    names(report) <- sprintf("log(%s)", what)
  }

  last <- cumsum(len)
  first <- c(0L, last[-length(last)]) + 1L
  x <- rep.int(NA_real_, length(time))
  ix <- list(
    interval = -first,
    cumulative = seq_along(time),
    rt = seq_along(time)
  )

  res <- data.frame(
    var = gl(length(report), length(time), labels = names(report)),
    ts = rep.int(frame_windows$ts[subset], len),
    window = rep.int(frame_windows$window[subset], len),
    time = time,
    estimate = unlist1(Map(replace, list(x), ix[what], report))
  )
  if (se) {
    res$se <- unlist(Map(replace, list(x), ix[what], report_se))
  }
  if (!log) {
    res$estimate <- exp(res$estimate)
    levels(res$var) <- what
  }
  attr(res, "se") <- se
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
#'   A logical flag. If \code{FALSE}, then confidence intervals
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
#' the \link{levels} of variable \code{var} modified accordingly.
#'
#' \code{level} is retained as an attribute.
#'
#' @examples
#' object <- egf_cache("predict-egf-1.rds")
#' confint(object, log = TRUE)
#' confint(object, log = FALSE)
#'
#' @export
confint.egf_predict <- function(object, parm, level = 0.95, log = TRUE, ...) {
  if (!isTRUE(attr(object, "se"))) {
    stop(wrap(
      "'object' must supply log scale predicted values ",
      "and corresponding standard errors. ",
      "Retry with 'object = predict(<\"egf\" object>, log = TRUE, se = TRUE)'."
    ))
  }
  stopifnot(
    is_number_in_interval(level, 0, 1, "()"),
    is_true_or_false(log)
  )

  res <- data.frame(
    object[-match("se", names(object), 0L)],
    wald(estimate = object$estimate, se = object$se, level = level)
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
