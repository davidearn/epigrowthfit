#' Compute confidence intervals on point estimates
#'
#' @description
#' A method for obtaining confidence intervals on point
#' estimates of model parameters from "egf" objects.
#'
#' @param object
#'   An "egf" object.
#' @param parm
#'   A subset of `c("R0", "doubling_time", names(object$theta_fit))`
#'   listing parameters. `NULL` is equivalent to the entire set,
#'   with `"R0"` excluded if either of `breaks` and `probs` is `NULL`.
#' @param level
#'   A number in the interval \[0,1\] indicating a confidence level.
#' @param method
#'   One of `"linear"`, `"uniroot"`, and `"wald"` indicating
#'   a method for computing the confidence interval.
#' @param breaks,probs
#'   Arguments to [compute_R0()] specifying a binned generation interval
#'   distribution. Necessary if `"R0" %in% parm` and unused otherwise.
#' @param ...
#'   Optional arguments to [TMB::tmbprofile()] or [TMB::tmbroot()],
#'   used if `method = "linear"` or `method = "uniroot"`, respectively.
#'
#' @return
#' If `length(parm) = 1`, then a numeric vector of the form
#' `c(estimate, lower, upper)`.
#' If `length(parm) > 1`, then a data frame of the form
#' `data.frame(estimate, lower, upper, row.names = parm)`.
#'
#' @seealso [egf()], [TMB::tmbprofile()], [TMB::tmbroot()]
#' @name confint.egf
NULL

#' @rdname confint.egf
#' @export
#' @import stats
#' @import TMB
confint.egf <- function(object,
                        parm = NULL,
                        level = 0.95,
                        method = "linear",
                        breaks = NULL,
                        probs = NULL,
                        ...) {
  opt <- c("doubling_time", names(object$theta_fit))
  if (!(is.null(breaks) || is.null(probs))) {
    opt <- c("R0", opt)
  }
  if (is.null(parm)) {
    parm <- opt
  } else {
    check(parm,
      what = "character",
      len = c(1, Inf),
      "`parm` must be `NULL` or a nonempty character vector."
    )
    check(parm,
      opt = opt,
      "Elements of `parm` must be chosen from:\n",
      paste(sprintf("\"%s\"", opt), collapse = ", ")
    )
    parm <- unique(parm)
  }
  check(level,
    what = "numeric",
    len = 1,
    val = c(0, 1),
    "`level` must be a number in the interval [0,1]."
  )
  check(method,
    what = "character",
    len = 1,
    opt = c("linear", "uniroot", "wald"),
    "`method` must be one of \"linear\", \"uniroot\", \"wald\"."
  )

  if (length(parm) > 1) {
    f <- function(x) {
      confint(object, parm = x, level = level, method = method, ...)
    }
    return(t(mapply(f, x = unname(parm))))
  }
  if (parm == "doubling_time") {
    ci <- confint(object, parm = "r", level = level, method = method, ...)
    return(setNames(compute_doubling_time(ci[c(1, 3, 2)]), names(ci)))
  }
  if (parm == "R0") {
    ci <- confint(object, parm = "r", level = level, method = method, ...)
    return(compute_R0(ci, breaks, probs))
  }

  sdr <- summary(sdreport(object$madf_out), select = "fixed")
  log_parm <- paste0("log_", parm)
  if (method == "linear") {
    p <- tmbprofile(object$madf_out, name = log_parm, ...)
    log_lu <- as.numeric(confint(p, level = level))
  } else if (method == "uniroot") {
    log_lu <- unname(
      tmbroot(object$madf_out,
        name = log_parm,
        target = 0.5 * qchisq(level, 1),
        ...
      )
    )
  } else if (method == "wald") {
    q <- qchisq(level, df = 1)
    ese <- sdr[log_parm, c("Estimate", "Std. Error")]
    log_lu <- ese[1] + c(-1, 1) * sqrt(q) * ese[2]
  }

  estimate <- object$theta_fit[[parm]]
  setNames(c(estimate, exp(log_lu)), c("estimate", "lower", "upper"))
}


