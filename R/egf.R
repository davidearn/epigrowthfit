#' \loadmathjax
#' Fit a model of epidemic growth
#'
#' @description
#' Minimizes the negative log likelihood of an epidemic growth model
#' using [`nlminb()`][stats::nlminb()], [`nlm()`][stats::nlm()], or
#' one of the routines provided through [`optim()`][stats::optim()].
#' The negative log likelihood function is written as a C++ template,
#' and its gradient with respect to log-transformed parameters is
#' defined by automatic differentiation using package \pkg{TMB}.
#'
#' @param init An "egf_init" object specifying epidemic data,
#'   a fitting window, and initial parameter estimates.
#'   See [egf_init()].
#' @param method One of `"nlminb"`, `"nlm"`, `"Nelder-Mead"`,
#'   `"BFGS"`, `"L-BFGS-S"`, and `"CG"`,
#'   indicating an optimization algorithm.
#' @param ... Additional arguments to [`nlminb()`][stats::nlminb()],
#'   [`nlm()`][stats::nlm()], or [`optim()`][stats::optim()].
#'
#' @return
#' An "egf" object. A list containing copies of arguments
#' `init` and `method`, with these additional elements:
#'
#' \describe{
#'   \item{`theta_hat`}{A named numeric vector giving the optimizer's
#'     approximation of the maximum likelihood parameter vector.
#'   }
#'   \item{`log_theta_hat`}{Log-transformed `theta_hat`. Identical
#'     to `log(theta_hat)` but with `"log_"` prepended to the names.
#'   }
#'   \item{`nll`}{A numeric scalar giving the negative log likelihood
#'     of `log_theta_hat`.
#'   }
#'   \item{`nll_func`}{A closure with a numeric argument `log_theta`
#'     specifying a log-transformed parameter vector (default is
#'     `log_theta_hat`). Returns the negative log likelihood function
#'     evaluated at `log_theta`.
#'   }
#'   \item{`nll_grad`}{A closure with a numeric argument `log_theta`
#'     specifying a log-transformed parameter vector (default is
#'     `log_theta_hat`). Returns the gradient of the negative log
#'     likelihood with respect to log-transformed parameters
#'     evaluated at `log_theta`.
#'   }
#'   \item{`eval_cum_inc`}{A closure with numeric arguments `time`
#'     and `theta` (default is `theta_hat`), evaluating expected
#'     cumulative incidence at `time` days using parameter
#'     vector `theta`. Elements must be named as in `theta_hat`.
#'   }
#'   \item{`madf_out`}{The list output of [TMB::MakeADFun()].}
#'   \item{`optim_out`}{The list output of
#'     [`nlminb()`][stats::nlminb()],
#'     [`nlm()`][stats::nlm()], or
#'     [`optim()`][stats::optim()],
#'     depending on `method`.
#'   }
#'   \item{`call`}{The call to `egf()`, allowing the output to
#'     be updated using [`update()`][stats::update()].
#'   }
#' }
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
#' sim <- simulate(x, nsim = 5, time = time_obs)
#' plot(sim, inc = "interval")
#' plot(sim, inc = "cumulative")
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
  if (!inherits(init, "egf_init")) {
    stop("`init` must be an \"egf_init\" object.")
  }
  valid_methods <- c("nlminb", "nlm", "Nelder-Mead", "BFGS", "L-BFGS-S", "CG")
  if (!is.character(method) || length(method) != 1 || !method %in% valid_methods) {
    warning("Invalid `method`, using `\"nlminb\"` instead.", call. = FALSE)
    method <- "nlminb"
  }

  ## Construct a call to `MakeADFun()`
  first <- init$first
  last <- init$last
  madf_data <- list(
    t = init$time[first:(last+1)] - init$time[first],
    x = init$cases[first:last],
    curve_flag = match(init$curve, c("exponential", "logistic", "richards")) - 1,
    baseline_flag = 1 * init$include_baseline,
    distr_flag = match(init$distr, c("pois", "nbinom")) - 1,
    method = method
  )
  log_par <- paste0("log_", c("r", "c0", "K", "thalf", "p", "nbdisp", "b"))
  madf_parameters <- as.list(init$log_theta0)
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
  ## super assignment in definition
  nll_func <- function(log_theta = log_theta_hat) {
    madf_out$fn(log_theta)
  }
  nll_grad <- function(log_theta = log_theta_hat) {
    madf_out$gr(log_theta)
  }

  ## Define a closure that evaluates expected
  ## cumulative incidence at desired time points
  eval_cum_inc <- function(time, theta = theta_hat) {
    ## Baseline for wave
    c1 <- if (first > 1) sum(init$cases[1:(first-1)]) else 0
    ## Excess following growth model
    c2 <- eval_model(time - init$time[first],
      curve = init$curve,
      include_baseline = init$include_baseline,
      theta = theta
    )
    c1 + c2
  }


  out <- list(
    init = init,
    method = method,
    theta_hat = theta_hat,
    log_theta_hat = log_theta_hat,
    nll = nll,
    nll_func = nll_func,
    nll_grad = nll_grad,
    eval_cum_inc = eval_cum_inc,
    madf_out = madf_out,
    optim_out = optim_out,
    call = match.call()
  )
  structure(out, class = c("egf", "list"))
}
