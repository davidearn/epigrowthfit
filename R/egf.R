#' \loadmathjax
#' Fit a model of epidemic growth
#'
#' @description
#' Minimizes the negative log likelihood of an epidemic growth model
#' using [stats::nlminb()], [stats::nlm()], of one of the routines
#' provided through [stats::optim()]. The negative log likelihood
#' function is written as a C++ template, and its gradient with
#' respect to log-transformed parameters is defined by automatic
#' differentiation using package \pkg{TMB}.
#'
#' @param init An "egf_init" object specifying epidemic data,
#'   a fitting window, and initial parameter estimates.
#'   See [egf_init()].
#' @param method One of `"nlminb"`, `"nlm"`, `"Nelder-Mead"`, `"BFGS"`,
#'   `"L-BFGS-S"`, and `"CG"`, indicating an optimization algorithm.
#' @param nbdisp_tol A positive number defining a threshold on the
#'   fitted value of the negative binomial dispersion parameter.
#'   Used only if `init$distr = "nbinom"`. See Details.
#' @param ... Additional arguments to [stats::nlminb()], [stats::nlm()],
#'   or [stats::optim()].
#'
#' @return
#' An "egf" object. A list containing copies of arguments
#' `init` and `method`, with these additional elements:
#'
#' \describe{
#'   \item{`theta_fit`}{A named numeric vector giving the optimizer's
#'     approximation of the maximum likelihood parameter vector.
#'   }
#'   \item{`log_theta_fit`}{Log-transformed `theta_fit`. Identical
#'     to `log(theta_fit)` but with `"log_"` prepended to the names.
#'   }
#'   \item{`nll`}{A numeric scalar giving the negative log likelihood
#'     of `log_theta_fit`.
#'   }
#'   \item{`nll_func`}{A closure with a numeric argument `log_theta`
#'     specifying a log-transformed parameter vector (default is
#'     `log_theta_fit`). Returns the negative log likelihood function
#'     evaluated at `log_theta`.
#'   }
#'   \item{`nll_grad`}{A closure with a numeric argument `log_theta`
#'     specifying a log-transformed parameter vector (default is
#'     `log_theta_fit`). Returns the gradient of the negative log
#'     likelihood with respect to log-transformed parameters evaluated
#'     at `log_theta`.
#'   }
#'   \item{`eval_cum_inc`}{A closure with numeric arguments `time`
#'     and `theta` (default is `theta_fit`), evaluating expected
#'     cumulative incidence at `time` days using parameter
#'     vector `theta`. Elements must be named as in `theta_fit`.
#'   }
#'   \item{`madf_out`}{The list output of [TMB::MakeADFun()].}
#'   \item{`optim_out`}{The list output of [stats::nlminb()],
#'     [stats::nlm()], or [stats::optim()], depending on `method`.
#'   }
#'   \item{`large_nbdisp_flag`}{A logical scalar. If `TRUE`, then
#'     the fitted value of the negative binomial dispersion parameter
#'     (i.e., `theta_fit[["nbdisp"]]`) exceeds the threshold defined
#'     by `nbdisp_tol`. Omitted if `init$distr != "nbinom"`.
#'     See Details.
#'   }
#'   \item{`call`}{The call to `egf()`, allowing the output to
#'     be updated using [stats::update()].
#'   }
#' }
#'
#' @details
#' If `theta_fit[["nbdisp"]]` exceeds
#'
#' `nbdisp_tol * max(diff(eval_cum_inc(init$time[init$first:(init$last+1)])))`
#'
#' then `egf()` will issue a warning suggesting to refit
#' using a Poisson model by running
#' `update(object, init = update(init, distr = "pois"))`,
#' where `object` is the "egf" object returned by `egf()`.
#' This behaviour is discussed in the package vignette,
#' accessible with `vignette("epigrowthfit-vignette")`.
#'
#' @examples
#' data(canadacovid)
#' ontario <- na.omit(subset(canadacovid, province == "ON"))
#' init <- egf_init(
#'   date = ontario$date,
#'   cases = ontario$new_confirmations[-1],
#'   curve = "logistic",
#'   distr = "nbinom"
#' )
#' x <- egf(init)
#' print(x)
#' coef(x, log = FALSE)
#' coef(x, log = TRUE)
#' plot(x, inc = "interval")
#' plot(x, inc = "cumulative")
#' time_obs <- init$time
#' time_pred <- seq(min(time_obs), max(time_obs), by = median(diff(time_obs)))
#' pred <- predict(x, time = time_pred)
#' sim <- simulate(x, nsim = 6, time = time_obs)
#' plot(sim, inc = "interval")
#' plot(sim, inc = "cumulative")
#'
#' @references
#' \insertRef{Ma+14}{epigrowthfit}
#'
#' \insertRef{Earn+20}{epigrowthfit}
#'
#' @seealso [egf_init()], [confint.egf()], [coef.egf()], [print.egf()],
#'   [predict.egf()], [simulate.egf()], [plot.egf()]
#' @export
#' @import stats
#' @importFrom TMB MakeADFun
#' @useDynLib epigrowthfit
egf <- function(init, method = "nlminb", nbdisp_tol = 100, ...) {
  if (!inherits(init, "egf_init")) {
    stop("`init` must be an \"egf_init\" object.")
  }
  m <- c("nlminb", "nlm", "Nelder-Mead", "BFGS", "L-BFGS-S", "CG")
  if (!is.character(method) || length(method) != 1 || !method %in% m) {
    warning("Invalid `method`, using `\"nlminb\"` instead.", call. = FALSE)
    method <- "nlminb"
  }
  if (init$distr == "nbinom") {
    if (!is.numeric(nbdisp_tol) || length(nbdisp_tol) != 1 ||
        !isTRUE(nbdisp_tol > 0)) {
      stop("`nbdisp_tol` must be a positive number.")
    }
  }

  ## Construct a call to `MakeADFun()`
  madf_data <- list(
    t = init$time[init$first:(init$last+1)] - init$time[init$first],
    x = init$cases[init$first:init$last],
    curve_flag = match(init$curve, c("exponential", "logistic", "richards")) - 1,
    baseline_flag = 1 * init$include_baseline,
    distr_flag = match(init$distr, c("pois", "nbinom")) - 1,
    method = method
  )
  log_par <- paste0("log_", c("r", "c0", "K", "thalf", "p", "nbdisp", "b"))
  madf_parameters <- as.list(init$log_theta_init)
  log_par_unused <- setdiff(log_par, names(madf_parameters))
  madf_parameters[log_par_unused] <- NA_real_
  madf_map <- setNames(
    replicate(length(log_par_unused), factor(NA), simplify = FALSE),
    log_par_unused
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

  ## Define a closure that evaluates expected
  ## cumulative incidence at desired time points
  eval_cum_inc <- function(time, theta = theta_fit) {
    ## Baseline for wave
    c1 <- if (init$first > 1) sum(init$cases[1:(init$first-1)]) else 0
    ## Excess following growth model
    c2 <- eval_model(time - init$time[init$first],
      curve = init$curve,
      include_baseline = init$include_baseline,
      theta = theta
    )
    c1 + c2
  }

  ## Warn if `nbdisp` exceeds threshold
  large_nbdisp_flag <- FALSE
  if (init$distr == "nbinom") {
    nbdisp_threshold <- max(diff(eval_cum_inc(init$time[init$first:(init$last+1)])))
    if (theta_fit[["nbdisp"]] > nbdisp_threshold) {
      warning("`nbdisp` exceeds threshold, refit with Poisson by running:\n\n",
              "update(object, init = update(init, distr = \"pois\"))",
              call. = FALSE)
      large_nbdisp_flag <- !large_nbdisp_flag
    }
  }


  out <- list(
    init = init,
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
  if (init$distr != "nbinom") {
    out$large_nbdisp_flag <- NULL
  }
  structure(out, class = c("egf", "list"))
}
