#' \loadmathjax
#' Fit a model of epidemic growth
#'
#' @description
#' Minimizes the negative log likelihood of an epidemic growth model
#' using [stats::nlminb()], [stats::nlm()], or one of the routines
#' provided through [stats::optim()]. The negative log likelihood as a
#' function of log-transformed parameters is written as a C++ template,
#' and its derivatives with respect to log-transformed parameters are
#' computed by automatic differentiation using package \pkg{TMB}.
#'
#' @param init
#'   An "egf_init" object specifying epidemic data, a fitting window,
#'   and initial parameter estimates. See [egf_init()].
#' @param method
#'   One of `"nlminb"`, `"nlm"`, `"Nelder-Mead"`, `"BFGS"`,
#'   `"L-BFGS-S"`, and `"CG"`, indicating an optimization algorithm.
#' @param na_action
#'   One of `"fail"` and `"exclude"`, indicating how missing values
#'   in the fitting window are handled.
#' @param nbdisp_tol
#'   A positive number defining a threshold on the fitted value
#'   of the negative binomial dispersion parameter. Used only if
#'   `init$distr = "nbinom"`.
#' @param ...
#'   Additional arguments to [stats::nlminb()], [stats::nlm()],
#'   or [stats::optim()], depending on `method`.
#'
#' @return
#' An "egf" object. A list containing copies of `init` elements
#' `data`, `window`, `first`, `last`, `curve`, `distr`, and
#' `include_baseline`; a copy of argument `method`; and these
#' additional elements:
#'
#' \describe{
#'   \item{`theta_fit`}{
#'     A named numeric vector giving the optimizer's approximation
#'     of the maximum likelihood parameter vector.
#'   }
#'   \item{`log_theta_fit`}{
#'     Log-transformed `theta_fit`. Identical to `log(theta_fit)`
#'     but with `"log_"` prepended to the names.
#'   }
#'   \item{`nll`}{
#'     A numeric scalar giving the negative log likelihood of
#'     `log_theta_fit`.
#'   }
#'   \item{`nll_func`}{
#'     A closure with a numeric argument `log_theta` specifying a
#'     log-transformed parameter vector (default is `log_theta_fit`).
#'     Returns the negative log likelihood function evaluated at
#'     `log_theta`.
#'   }
#'   \item{`nll_grad`}{
#'     A closure with a numeric argument `log_theta` specifying a
#'     log-transformed parameter vector (default is `log_theta_fit`).
#'     Returns the gradient of the negative log likelihood with
#'     respect to log-transformed parameters evaluated at `log_theta`.
#'   }
#'   \item{`eval_cum_inc`}{
#'     A closure with numeric arguments `time` and `theta`
#'     (default is `theta_fit`) evaluating expected cumulative
#'     incidence at `time` days since `data$date[first]` conditional
#'     on parameter vector `theta`. Elements of `theta` must be
#'     named as in `theta_fit`. `predict(object, time)` should
#'     be used instead of `object$eval_cum_inc(time)`.
#'   }
#'   \item{`madf_out`}{The list output of [TMB::MakeADFun()].}
#'   \item{`optim_out`}{
#'     The list output of [stats::nlminb()], [stats::nlm()],
#'     or [stats::optim()], depending on `method`.
#'   }
#'   \item{`large_nbdisp_flag`}{
#'     A logical scalar. If `TRUE`, then the fitted value of
#'     the negative binomial dispersion parameter
#'     (i.e., `theta_fit[["nbdisp"]]`) exceeds the threshold
#'     defined by `nbdisp_tol`. Omitted if `distr != "nbinom"`.
#'   }
#'   \item{`call`}{
#'     The call to `egf()`, allowing the output to be updated
#'     using [stats::update()].
#'   }
#' }
#'
#' @details
#' The data used in likelihood calculation are specified by
#' `x = with(data, cases[(first+1):last])`.
#' If `x` has missing values and `na_action = "fail"`,
#' then `egf()` will stop with an error.
#' If `x` has missing values and `na_action = "exclude"`,
#' then the likelihood function will be defined as a product
#' over indices without missing values, unless `sum(!is.na(x))`
#' is less than the number of model parameters, in which case
#' `egf()` will stop with an error.
#'
#' If `distr = "nbinom"` and `theta_fit[["nbdisp"]]`
#' exceeds `nbdisp_tol * max(diff(eval_cum_inc(time)))`,
#' where `time = data$time[first:last]`, then `egf()`
#' will issue a warning suggesting to refit using a
#' Poisson model by running
#' `update(object, init = update(init, distr = "pois"))`.
#' This behaviour is discussed in the package vignette,
#' accessible with `vignette("epigrowthfit-vignette")`.
#'
#' @examples
#' data(canadacovid)
#' ontario <- subset(canadacovid, province == "ON")
#' x1 <- egf_init(new_confirmed ~ date,
#'   data = ontario,
#'   curve = "logistic",
#'   distr = "nbinom",
#'   first = 27,
#'   last = 77,
#'   peak = 77
#' )
#' x2 <- update(x1, first = 211, last = 236, peak = "2020-10-10")
#' y1 <- egf(x1)
#' y2 <- egf(x2)
#' print(y1)
#' coef(y1, log = FALSE)
#' coef(y1, log = TRUE)
#' confint(y1, parm = "doubling_time", level = 0.95, method = "linear")
#' predict(y1)
#' s <- simulate(y1, nsim = 6)
#' plot(s, inc = "interval")
#' plot(s, inc = "cumulative")
#' plot(y1, inc = "interval")
#' plot(y2, inc = "interval", add = TRUE)
#' plot(y1, inc = "cumulative")
#' plot(y2, inc = "cumulative", add = TRUE)
#'
#' @references
#' \insertRef{Ma+14}{epigrowthfit}
#'
#' \insertRef{Earn+20}{epigrowthfit}
#'
#' @seealso [egf_init()], [confint.egf()], [plot.egf()]
#' @export
#' @import stats
#' @importFrom TMB MakeADFun
#' @useDynLib epigrowthfit
egf <- function(init, method = "nlminb", na_action = "exclude",
                nbdisp_tol = 100, ...) {
  check(init,
    what = "egf_init",
    "`init` must be an \"egf_init\" object."
  )
  ok <- check(method,
    what = "character",
    len = 1,
    opt = c("nlminb", "nlm", "Nelder-Mead", "BFGS", "L-BFGS-S", "CG"),
    action = "warn",
    "Invalid `method`, using \"nlminb\" instead."
  )
  if (!ok) {
    method <- "nlminb"
  }
  ok <- check(na_action,
    what = "character",
    len = 1,
    opt = c("fail", "exclude"),
    action = "warn",
    "Invalid `na_action`, using \"exclude\" instead."
  )
  if (!ok) {
    na_action <- "exclude"
  }
  if (init$distr == "nbinom") {
    check(nbdisp_tol,
      what = "numeric",
      len = 1,
      val = c(0, Inf),
      no = is.na,
      "`nbdisp_tol` must be a non-negative number."
    )
  }

  ## Construct `MakeADFun()` input
  madf_data <- with(init,
    list(
      t = data$time[first:last],
      x = data$cases[(first+1):last],
      curve_flag = match(curve, c("exponential", "logistic", "richards")) - 1,
      baseline_flag = 1 * include_baseline,
      distr_flag = match(distr, c("pois", "nbinom")) - 1,
      method = method
    )
  )
  log_par <- paste0("log_", c("r", "c0", "K", "thalf", "p", "nbdisp", "b"))
  madf_parameters <- as.list(init$log_theta_init)
  log_par_unused <- setdiff(log_par, names(madf_parameters))
  madf_parameters[log_par_unused] <- NA_real_
  madf_map <- setNames(
    replicate(length(log_par_unused), factor(NA), simplify = FALSE),
    log_par_unused
  )

  ## Perform checks for `NA`
  npar <- with(init, length(theta_init))
  if (na_action == "fail") {
    check(madf_data$x,
      no = anyNA,
      "Fitting window contains missing values."
    )
  } else if (na_action == "exclude") {
    check(madf_data$x,
      no = function(x) sum(!is.na(x)) < npar,
      "Fitting window contains insufficient data."
    )
  }

  ## Call `MakeADFun()` to define the negative log likelihood
  ## function and its gradient
  madf_out <- MakeADFun(
    data = madf_data,
    parameters = madf_parameters,
    map = madf_map,
    DLL = "epigrowthfit",
    silent = TRUE
  )

  ## Call `nlminb()`, `nlm()`, or `optim()` (depending on `method`)
  ## to minimize the negative log likelihood function
  if (method == "nlminb") {
    optim_out <- nlminb(
      start = madf_out$par,
      objective = madf_out$fn,
      gradient = madf_out$gr,
      ...
    )
    log_theta_fit <- optim_out$par
    nll <- optim_out$objective
  } else if (method == "nlm") {
    optim_out <- nlm(
      f = structure(madf_out$fn, gradient = madf_out$gr),
      p = madf_out$par,
      ...
    )
    log_theta_fit <- optim_out$estimate
    nll <- optim_out$minimum
  } else {
    optim_out <- optim(
      par = madf_out$par,
      fn = madf_out$fn,
      gr = madf_out$gr,
      method = method,
      ...
    )
    log_theta_fit <- optim_out$par
    nll <- optim_out$value
  }

  ## Inverse-link transform the fitted parameter vector
  theta_fit <- setNames(
    exp(log_theta_fit),
    sub("log_", "", names(log_theta_fit))
  )

  ## Define wrappers for the negative log likelihood
  ## function and its gradient to circumvent unexpected
  ## super assignment in definition
  nll_func <- function(log_theta = log_theta_fit) {
    madf_out$fn(log_theta)
  }
  nll_grad <- function(log_theta = log_theta_fit) {
    madf_out$gr(log_theta)
  }

  ## Define a closure that evaluates expected cumulative incidence
  ## at times in days since the start of the fitting window
  eval_cum_inc <- function(time, theta = theta_fit) {
    eval_model(time,
      curve = init$curve,
      include_baseline = init$include_baseline,
      theta = theta
    )
  }

  ## Warn if `nbdisp` exceeds threshold
  large_nbdisp_flag <- FALSE
  if (init$distr == "nbinom") {
    nbdisp_threshold <- max(diff(eval_cum_inc(madf_data$t)))
    if (theta_fit[["nbdisp"]] > nbdisp_threshold) {
      warning("`nbdisp` exceeds threshold, refit with Poisson by running:\n\n",
              "update(object, init = update(init, distr = \"pois\"))",
              call. = FALSE)
      large_nbdisp_flag <- !large_nbdisp_flag
    }
  }


  out1 <- init[c("data", "window", "first", "last",
                 "curve", "distr", "include_baseline")]
  out2 <- list(
    method = method,
    theta_fit = theta_fit,
    log_theta_fit = log_theta_fit,
    nll = nll,
    nll_func = nll_func,
    nll_grad = nll_grad,
    eval_cum_inc = eval_cum_inc,
    madf_out = madf_out,
    optim_out = optim_out,
    large_nbdisp_flag = large_nbdisp_flag,
    call = match.call()
  )
  if (out1$distr != "nbinom") {
    out2$large_nbdisp_flag <- NULL
  }
  structure(c(out1, out2), class = c("egf", "list"))
}
