#' \loadmathjax
#' Fit models of epidemic growth
#'
#' @description
#' A description of function `egf()`.
#'
#' @param x An "egf_init" object, specifying epidemic data,
#'   a fitting window, and initial parameter estimates.
#'   See [egf_init()].
#' @param method One of `"Nelder-Mead"`, `"BFGS"`, `"CG"`,
#'   `"L-BFGS-S"`, and `"SANN"`, indicating a method used
#'   to minimize the negative log likelihood function.
#'   Passed directly to [stats::optim()].
#' @param ... Additional arguments to [stats::optim()].
#'
#' @return
#' An "egf" object.
#'
#' @details
#' ## Phenomenological model
#' Let \mjseqn{x(t)} be the expected number of cases observed
#' up to time \mjseqn{t}, i.e., expected cumulative incidence.
#' \mjseqn{x(t)} is modeled by an exponential, logistic, or
#' Richards curve, plus a linear term representing baseline
#' mortality, relevant if \mjseqn{x(t)} counts all causes deaths
#' instead of disease deaths (and otherwise excluded).
#'
#' ### Exponential
#'
#' \mjsdeqn{x(t) = b t + x_0 e^{r t}}
#'
#' ### Logistic
#'
#' \mjsdeqn{x(t) = b t + \frac{K}{1 + \big(\frac{K}{x(0)} - 1\big) e^{-r t}}}
#'
#' reparametrized as
#'
#' \mjsdeqn{x(t) = b t + \frac{K}{1 + e^{-r (t - t_\text{half}}}}
#'
#' where \mjseqn{t_\text{half}} satisfies
#' \mjseqn{x(t_\text{half}) - b t_\text{half} = \frac{K}{2}}.
#'
#' ### Richards
#'
#' \mjsdeqn{x(t) = b t + \frac{K}{\big\lbrack 1 + \big(\big(\frac{K}{x(0)}\big)^p - 1\big) e^{-r p t} \big\rbrack^{1/p}}}
#'
#' reparametrized as
#'
#' \mjsdeqn{x(t) = b t + \frac{K}{\big\lbrack 1 + (2^p - 1) * e^{-r p (t - t_\text{half}} \big\rbrack}^{1/p}}
#'
#' where \mjseqn{t_\text{half}} satisfies
#' \mjseqn{x(t_\text{half}) - b t_\text{half} = \frac{K}{2}}.
#'
#' ## Observation model
#' Let \mjseqn{Y(t_1,t_2)} be the number of cases observed between times
#' \mjseqn{t_1} and \mjseqn{t_2 > t_1}, i.e., observed interval incidence.
#' \mjseqn{Y(t_1,t_2)} is modeled as either a Poisson-distributed random
#' variable with mean \mjseqn{x(t_2) - x(t_1)},
#'
#' \mjsdeqn{Y(t_1,t_2) \sim \mathrm{Poisson}\big(x(t_2) - x(t_1)\big)\,,}
#'
#' or as negative binomial-distributed random variable with mean
#' \mjseqn{x(t_2) - x(t_1)} and dispersion \mjseqn{k},
#'
#' \mjsdeqn{Y(t_1,t_2) \sim \mathrm{NegativeBinomial}\big(x(t_2) - x(t_1), k\big)\,.}
#'
#' @examples
#' data(canadacovid)
#' ontario <- na.omit(subset(canadacovid, province == "ON"))
#' init <- egf_init(
#'   date = ontario$date,
#'   cases = ontario$new_confirmations[-1],
#'   curve = "richards",
#'   distr = "nbinom"
#' x <- egf(init)
#' plot(x, inc = "interval")
#' plot(x, inc = "cumulative")
#'
#' @export
#' @importFrom TMB MakeADFun
#' @importFrom stats setNames optim
#' @useDynLib epigrowthfit
egf <- function(x, method = "BFGS", ...) {
  if (missing(x)) {
    stop("Missing argument `x`.")
  } else if (!inherits(x, "egf_init")) {
    stop("`x` must be an \"egf_init\" object.")
  }
  if (!is.character(method) || length(method) != 1 ||
        !method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-S", "SANN")) {
    stop("`method` must be an element of\n",
         "`c(\"Nelder-Mead\", \"BFGS\", \"CG\", \"L-BFGS-S\", \"SANN\")`.")
  }

  f <- x$first
  l <- x$last
  madf_data <- list(
    t = x$time[f:(l+1)],
    x = x$cases[f:l],
    curve_flag = match(x$curve, c("exponential", "logistic", "richards")) - 1,
    baseline_flag = 1 * x$include_baseline,
    distr_flag = match(x$distr, c("pois", "nbinom")) - 1
  )

  theta0 <- x$theta0
  pars_unused <- setdiff(
    c("r", "x0", "K", "thalf", "p", "b", "nbdisp"),
    names(theta0)
  )
  theta0[pars_unused] <- NA_real_
  madf_parameters <- setNames(
    lapply(theta0, log),
    paste0("log_", names(theta0))
  )

  madf_map <- setNames(
    replicate(length(pars_unused), factor(NA), simplify = FALSE),
    paste0("log_", pars_unused)
  )

  madf_out <- MakeADFun(
    data = madf_data,
    parameters = madf_parameters,
    map = madf_map,
    DLL = "epigrowthfit",
    silent = TRUE
  )
  optim_out <- optim(
    par = madf_out$par,
    fn = madf_out$fn,
    gr = madf_out$gr,
    method = method,
    ...
  )

  theta_mle <- exp(optim_out$par)
  names(theta_mle) <- sub("log_", "", names(theta_mle))
  eval_cum_inc <- function(time, theta = theta_mle) {
    with(as.list(theta), {
      cum_inc <- switch(x$curve,
        exponential = x0 * exp(r * time),
        logistic = K / (1 + exp(-r * (time - thalf))),
        richards = K / (1 + (2^p - 1) * exp(-r * p * (time - thalf)))^(1 / p)
      )
      if (x$include_baseline) {
        cum_inc <- b * time + cum_inc
      }
      cum_inc
    })
  }
  eval_int_inc <- function(time, theta = theta_mle) {
    diff(eval_cum_inc(time, theta))
  }

  out <- list(
    theta_mle = theta_mle,
    eval_cum_inc = eval_cum_inc,
    eval_int_inc = eval_int_inc,
    init = x,
    madf_out = madf_out,
    optim_out = optim_out
  )
  structure(out, class = c("egf", "list"))
}
