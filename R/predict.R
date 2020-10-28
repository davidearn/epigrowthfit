##' Evaluate models of epidemic growth
#'
#' @description
#' Methods for evaluating incidence models specified
#' by objects of class "egf_init" or "egf".
#'
#' @param object An "egf_init" or "egf" object.
#' @param time A numeric vector listing increasing time points,
#'   expressed as numbers of days since `object$date[1]`
#'   for "egf_init" objects and since `object$init$date[1]`
#'   for "egf" objects. Missing values are not tolerated.
#' @param ... Unused optional arguments.
#'
#' @return
#' A list with numeric elements:
#'
#' \describe{
#'   \item{`time`}{Matches argument.}
#'   \item{`refdate`}{Matches `object$date[1]` or `object$init$date[1]`.}
#'   \item{`cum_inc`}{Expected cumulative incidence at time points
#'     `time`. Equal to `object$eval_cum_inc(time)`.
#'   }
#'   \item{`int_inc`}{Expected interval incidence given interval
#'     endpoints `time`. Equal to `diff(object$eval_cum_inc(time))`
#'     if `length(time) >= 2` and omitted otherwise.
#'   }
#' }

#'
#' @seealso [egf_init()], [egf()]
#' @name predict.egf
NULL

#' @rdname predict.egf
#' @export
predict.egf_init <- function(object, time = object$time, ...) {
  if (!is.numeric(time) || length(time) == 0) {
    stop("`time` must be numeric and have nonzero length.")
  } else if (anyNA(time)) {
    stop("`time` must not have missing values.")
  } else if (!all(diff(time) > 0)) {
    stop("`time` must be increasing.")
  } else if (any(time < object$time[object$first] |
                 time > object$time[object$last+1])) {
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

#' @rdname predict.egf
#' @export
predict.egf <- function(object, time = object$init$time, ...) {
  if (!is.numeric(time) || length(time) == 0) {
    stop("`time` must be numeric and have nonzero length.")
  } else if (anyNA(time)) {
    stop("`time` must not have missing values.")
  } else if (!all(diff(time) > 0)) {
    stop("`time` must be increasing.")
  } else if (any(time < object$init$time[object$init$first] |
                 time > object$init$time[object$init$last+1])) {
    warning("There are elements of `time` outside of the fitting window.",
            call. = FALSE)
  }


  out <- list(
    time = time,
    refdate = object$init$date[1],
    cum_inc = object$eval_cum_inc(time)
  )
  if (length(time) > 1) {
    out$int_inc = diff(out$cum_inc)
  }
  out
}
