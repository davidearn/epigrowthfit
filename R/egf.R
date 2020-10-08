#' \loadmathjax
#' Fit models of epidemic growth
#'
#' @description
#' Minimizes the negative log likelihood of an epidemic growth model
#' using [`nlminb()`][stats::nlminb()], [`nlm()`][stats::nlm()], or
#' one of the routines provided through [`optim()`][stats::optim()].
#' The negative log likelihood function is written as a C++ template.
#' Its gradient with respect to parameters is computed by automatic
#' differentiation using package \pkg{TMB}.
#'
#' @param init An "egf_init" object, specifying epidemic data,
#'   a fitting window, and initial parameter estimates.
#'   See [egf_init()].
#' @param method One of `"nlminb"`, `"nlm"`, `"Nelder-Mead"`,
#'   `"BFGS"`, `"L-BFGS-S"`, and `"CG"`, indicating the method
#'   used to minimize the negative log likelihood function.
#' @param ... Additional arguments to [`nlminb()`][stats::nlminb()],
#'   [`nlm()`][stats::nlm()], or [`optim()`][stats::optim()].
#'
#' @return
#' An "egf" object. A list with elements:
#'
#' \describe{
#'   \item{`theta_hat`}{A named numeric vector specifying the
#'     optimizer's estimate of the minimizer of the negative
#'     log likelihood function.
#'   }
#'   \item{`log_theta_hat`}{Log-transformed `theta_hat`. Identical
#'     to `log(theta_hat)` with `"log_"` prepended to the names.
#'   }
#'   \item{`nll`}{A numeric scalar giving the negative log likelihood
#'     corresponding to `log_theta_hat`.
#'   }
#'   \item{`nll_fn`}{A closure with numeric argument `log_theta`
#'     specifying a log-transformed parameter vector (default
#'     is `log_theta_hat`). Returns the negative log likelihood
#'     function evaluated at `x`.
#'   }
#'   \item{`nll_gr`}{A closure like `nll_fn`, but returning the
#'     gradient of the negative log likelihood function evaluated
#'     at `x`.
#'   }
#'   \item{`cum_inc`}{A closure with numeric arguments `time`
#'     and `theta` (default is `theta_hat`), evaluating expected
#'     cumulative incidence at `time` using parameter values
#'     `theta`.
#'   }
#'   \item{`int_inc`}{A closure with numeric arguments `time`
#'     and `theta` (default is `theta_hat`), evaluating expected
#'     interval incidence given `time` as interval endpoints.
#'     Returns `diff(cum_inc(time, theta))`, hence the value
#'     has length `length(time)-1`.
#'   }
#'   \item{`init`}{Matches argument.}
#'   \item{`madf_out`}{The list output of [TMB::MakeADFun()].}
#'   \item{`optim_out`}{The list output of [`nlminb()`][stats::nlminb()],
#'     [`nlm()`][stats::nlm()], or [`optim()`][stats::optim()],
#'     depending on `method`.
#'   }
#'   \item{`call`}{The call to `egf()`, making the output
#'     reproducible with `eval(call)`.
#'   }
#' }
#'
#' @examples
#' data(canadacovid)
#' ontario <- na.omit(subset(canadacovid, province == "ON"))
#' init <- egf_init(
#'   date = ontario$date,
#'   cases = ontario$new_confirmations[-1],
#'   curve = "richards",
#'   distr = "nbinom"
#' )
#' x <- egf(init)
#' print(x)
#' coef(x, log = FALSE)
#' coef(x, log = TRUE)
#' time_obs <- as.numeric(ontario$date - ontario$date[1])
#' time_pred <- seq(min(time_obs), max(time_obs), by = median(diff(time_obs)))
#' predict(x, time = time_pred)
#' simulate(x, nsim = 4, time = time_obs)
#' plot(x, inc = "interval")
#' plot(x, inc = "cumulative")
#'
#' @references
#' \insertRef{Ma+14}{epigrowthfit}
#'
#' \insertRef{Earn+20}{epigrowthfit}
#'
#' @seealso [egf_init()], [methods for class "egf"][egf-methods]
#'
#' @export
#' @import stats
#' @importFrom TMB MakeADFun
#' @useDynLib epigrowthfit
egf <- function(init, method = "nlminb", ...) {
  if (missing(init)) {
    stop("Missing argument `init`.")
  } else if (!inherits(init, "egf_init")) {
    stop("`init` must be an \"egf_init\" object.")
  }
  valid_methods <- c("nlminb", "nlm", "Nelder-Mead", "BFGS", "L-BFGS-S", "CG")
  if (!is.character(method) || length(method) != 1 ||
        !method %in% valid_methods) {
    warning("Invalid `method`, using `\"nlminb\"` instead.")
    method <- "nlminb"
  }

  ## Construct a call to `MakeADFun()`
  f <- init$first
  l <- init$last
  madf_data <- list(
    t = init$time[f:(l+1)],
    x = init$cases[f:l],
    curve_flag = match(init$curve, c("exponential", "logistic", "richards")) - 1,
    baseline_flag = 1 * init$include_baseline,
    distr_flag = match(init$distr, c("pois", "nbinom")) - 1,
    method = method
  )
  par <- paste0("log_", c("r", "x0", "K", "thalf", "p", "b", "nbdisp"))
  madf_parameters <- as.list(init$log_theta0)
  par_unused <- setdiff(par, names(madf_parameters))
  madf_parameters[par_unused] <- NA_real_
  madf_map <- setNames(
    replicate(length(par_unused), factor(NA), simplify = FALSE),
    par_unused
  )

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
    log_theta_hat <- optim_out$par
    nll <- optim_out$objective
  } else if (method == "nlm") {
    optim_out <- nlm(
      f = structure(madf_out$fn, gradient = madf_out$gr),
      p = madf_out$par,
      ...
    )
    log_theta_hat <- optim_out$estimate
    nll <- optim_out$minimum
  } else {
    optim_out <- optim(
      par = madf_out$par,
      fn = madf_out$fn,
      gr = madf_out$gr,
      method = method,
      ...
    )
    log_theta_hat <- optim_out$par
    nll <- optim_out$value
  }

  ## Inverse-link transform the fitted parameter values
  theta_hat <- setNames(
    exp(log_theta_hat),
    sub("log_", "", names(log_theta_hat))
  )

  ## Define wrappers for the negative log likelihood
  ## function and its gradient to circumvent unexpected
  ## global assignment in their definition
  nll_fn <- function(log_theta = log_theta_hat) {
    madf_out$fn(log_theta)
  }
  nll_gr <- function(log_theta = log_theta_hat) {
    madf_out$gr(log_theta)
  }

  ## Define convenience functions for evaluating
  ## cumulative and interval incidence
  cum_inc <- function(time, theta = theta_hat) {
    x <- eval_inc(time,
      curve = init$curve,
      include_baseline = init$include_baseline,
      theta = theta
    )
    x$cum_inc
  }
  int_inc <- function(time, theta = theta_hat) {
    x <- eval_inc(time,
      curve = init$curve,
      include_baseline = init$include_baseline,
      theta = theta
    )
    x$int_inc
  }


  out <- list(
    theta_hat = theta_hat,
    log_theta_hat = log_theta_hat,
    nll = nll,
    nll_fn = nll_fn,
    nll_gr = nll_gr,
    cum_inc = cum_inc,
    int_inc = int_inc,
    init = init,
    madf_out = madf_out,
    optim_out = optim_out,
    call = match.call()
  )
  structure(out, class = c("egf", "list"))
}
