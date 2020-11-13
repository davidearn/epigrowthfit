#' Confidence intervals on point estimates
#'
#' @description
#' A method for obtaining confidence intervals on point estimates
#' of model parameters specified by "egf" objects.
#'
#' @param object
#'   An "egf" object.
#' @param parm
#'   An element of `names(object$theta_fit)`, or otherwise
#'   `"doubling_time"`, indicating a parameter for which to
#'   compute a confidence interval.
#' @param level
#'   A number in the interval \[0,1\] indicating a confidence level.
#' @param method
#'   One of `"linear"`, `"uniroot"`, and `"wald"` indicating
#'   a method with which to compute the confidence interval.
#' @param ... Unused optional arguments.
#'
#' @return
#' A numeric vector of the form `c(lower = a, upper = b)`,
#' where `a` and `b` are the left and right endpoints of
#' the confidence interval specified by `parm`, `level`,
#' and `method`.
#'
#' @seealso [egf()], [TMB::tmbprofile()], [TMB::tmbroot()]
#' @name confint.egf
NULL

#' @rdname confint.egf
#' @export
#' @import stats
#' @import TMB
confint.egf <- function(object,
                        parm = "doubling_time",
                        level = 0.95,
                        method = "linear",
                        ...) {
  check(parm,
    what = "character",
    len = 1,
    opt = c("doubling_time", names(object$theta_fit)),
    "`parm` must be an element of `names(object$theta_fit)`\n",
    "or \"doubling_time\"."
  )
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

  if (parm == "doubling_time") {
    ci <- confint(object, parm = "r", level = level, method = method)
    ci <- setNames(rev(compute_doubling_time(ci)), names(ci))
    return(ci)
  }

  sdr <- summary(sdreport(object$madf_out), select = "fixed")
  log_parm <- paste0("log_", parm)
  if (method == "linear") {
    cat <- function(...) {}
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

  setNames(exp(log_lu), c("lower", "upper"))
}


