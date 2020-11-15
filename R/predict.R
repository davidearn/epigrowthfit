##' Evaluate models of epidemic growth
#'
#' @description
#' Methods for evaluating incidence models specified
#' by objects of class "egf_init" or "egf".
#'
#' @param object
#'   An "egf_init" or "egf" object.
#' @param time
#'   A numeric vector listing increasing time points.
#'   Times must be expressed as numbers of days since
#'   `with(object, data$date[first])`.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A data frame with numeric variables:
#'
#' \describe{
#'   \item{`time`}{Matches argument.}
#'   \item{`cum_inc`}{
#'     Expected cumulative incidence at time points `time`.
#'     Equal to `object$eval_cum_inc(time)`.
#'   }
#'   \item{`int_inc`}{
#'     Expected interval incidence given interval endpoints `time`.
#'     Equal to `c(NA, diff(cum_inc))`.
#'   }
#' }
#'
#' Attribute `refdate = with(object, data$date[first])` gives
#' the reference date from which times in `time` are measured.
#'
#' @details
#' A warning will be issued if any elements of `time` fall
#' outside of the fitting window, indicating that the model
#' is being extrapolated. To avoid warnings, make sure that
#' elements of `time` are in the closed interval
#' `with(object, data$time[c(first, last)])`.
#'
#' @seealso [egf_init()], [egf()]
#' @name predict.egf
NULL

#' @rdname predict.egf
#' @export
predict.egf_init <- function(object,
                             time = with(object, data$time[window]),
                             ...) {
  ## Reuse `predict.egf()` machinery
  predict.egf(object, time)
}

#' @rdname predict.egf
#' @export
predict.egf <- function(object,
                        time = with(object, data$time[window]),
                        ...) {
  check(time,
    what = "numeric",
    len = c(1, Inf),
    "`time` must be numeric and have nonzero length."
  )
  check(time,
    no = anyNA,
    "`time` must not have missing values."
  )
  check(time,
    yes = function(x) all(diff(x) > 0),
    "`time` must be increasing."
  )
  check(time,
    val = with(object, data$time[c(first, last)]),
    action = "warn",
    "There are elements of `time` outside of the fitting window."
  )

  cum_inc_pred <- object$eval_cum_inc(time)
  int_inc_pred <- diff(cum_inc_pred)
  out <- data.frame(
    time = time,
    cum_inc = cum_inc_pred,
    int_inc = c(NA, int_inc_pred)
  )
  structure(out, refdate = with(object, data$date[first]))
}
