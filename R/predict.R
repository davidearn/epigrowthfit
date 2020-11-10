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
#'   `with(object, date[first])` for "egf_init" objects
#'   and since
#'   `with(object$init, date[first])` for "egf" objects.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A list with numeric elements:
#'
#' \describe{
#'   \item{`time`}{Matches argument.}
#'   \item{`refdate`}{
#'     Matches `with(object, date[first])`
#'     or `with(object$init, date[first])`.
#'   }
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
#' @seealso [egf_init()], [egf()]
#' @name predict.egf
NULL

#' @rdname predict.egf
#' @export
predict.egf_init <- function(object,
                             time = with(object, time[first:last]),
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
    val = with(object, time[c(first, last)]),
    action = "warn",
    "There are elements of `time` outside of the fitting window."
  )

  cum_inc_pred <- object$eval_cum_inc(time)
  int_inc_pred <- diff(cum_inc_pred)
  list(
    time = time,
    refdate = with(object, date[first]),
    cum_inc = cum_inc_pred,
    int_inc = c(NA, int_inc_pred)
  )
}

#' @rdname predict.egf
#' @export
predict.egf <- function(object,
                        time = with(object$init, time[first:last]),
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
    val = with(object$init, time[c(first, last)]),
    action = "warn",
    "There are elements of `time` outside of the fitting window."
  )

  cum_inc_pred <- object$eval_cum_inc(time)
  int_inc_pred <- diff(cum_inc_pred)
  list(
    time = time,
    refdate = with(object$init, date[first]),
    cum_inc = cum_inc_pred,
    int_inc = c(NA, int_inc_pred)
  )
}
