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
#' @param se
#'   A logical scalar. If `TRUE`, then standard errors on
#'   predicted incidence are computed via the delta method.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' If `se = FALSE`, then a data frame with numeric variables:
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
#' If `se = TRUE`, then an "egf_pred" object ...
#'
#' Attribute `refdate = with(object, data$date[first])` gives
#' the reference date from which times in `time` are measured.
#'
#' @details
#' A warning will be issued if any elements of `time` fall
#' outside of the fitting window, i.e., if `predict()` is
#' asked to extrapolate the model. To avoid warnings, make
#' sure that elements of `time` are in the closed interval
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
  predict.egf(object, time, se = FALSE)
}

#' @rdname predict.egf
#' @export
#' @importFrom stats setNames
#' @importFrom TMB MakeADFun sdreport
predict.egf <- function(object,
                        time = with(object, data$time[window]),
                        se = FALSE,
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

  if (se) {

    ## FIXME: DRY!
    ## Up to `MakeADFun()` call, much is repeated from `egf()`...
    madf_data <- with(object,
      list(
        t = data$time[first:last],
        x = data$cases[(first+1):last],
        curve_flag = match(curve, c("exponential", "logistic", "richards")) - 1L,
        distr_flag = match(distr, c("pois", "nbinom")) - 1L,
        baseline_flag = 1L * include_baseline,
        predict_flag = 1L,
        t_new = time,
        x_new = rep(NA_real_, length(time) - 1)
      )
    )
    log_par <- paste0("log_", c("r", "c0", "K", "thalf", "p", "nbdisp", "b"))
    madf_parameters <- as.list(object$log_theta_fit)
    log_par_unused <- setdiff(log_par, names(madf_parameters))
    madf_parameters[log_par_unused] <- NA_real_
    madf_map <- setNames(
      replicate(length(log_par_unused), factor(NA), simplify = FALSE),
      log_par_unused
    )
    madf_out <- MakeADFun(
      data = madf_data,
      parameters = madf_parameters,
      map = madf_map,
      DLL = "epigrowthfit",
      silent = TRUE
    )

    madf_out$fn(object$log_theta_fit) # ?!
    sdr <- summary(sdreport(madf_out), select = "report")
    sdr_split <- split(as.data.frame(sdr), factor(rownames(sdr)))
    names(sdr_split) <- sub("_new", "", names(sdr_split), fixed = TRUE)
    sdr_split$log_int_inc <- rbind(NA, sdr_split$log_int_inc)
    out <- lapply(sdr_split, function(x) {
      names(x) <- c("estimate", "se")
      row.names(x) <- NULL
      cbind(time, x)
    })
    class(out) <- c("egf_pred", "list")

  } else {

    cum_inc_pred <- object$eval_cum_inc(time)
    int_inc_pred <- diff(cum_inc_pred)
    out <- data.frame(
      time = time,
      cum_inc = cum_inc_pred,
      int_inc = c(NA, int_inc_pred)
    )

  }

  structure(out, refdate = with(object, data$date[first]))
}
