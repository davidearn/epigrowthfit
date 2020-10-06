#' \loadmathjax
#' Fit models of epidemic growth
#'
#' @description
#' A description of function `egf()`.
#'
#' @param init An "egf_init" object, specifying epidemic data,
#'   a fitting window, and initial parameter estimates.
#'   See [egf_init()].
#' @param method One of `"Nelder-Mead"`, `"BFGS"`, `"CG"`,
#'   `"L-BFGS-S"`, and `"SANN"`, indicating a method used
#'   to minimize the negative log likelihood function.
#'   Passed directly to [stats::optim()].
#' @param ... Additional arguments to [stats::optim()].
#'
#' @return
#' An "egf" object. A list with elements
#'
#' \describe{
#'   \item{`theta_mle`}{A list containing the maximum likelihood
#'     estimate of the model parameters.
#'   }
#'   \item{`eval_cum_inc`}{A closure with arguments `time` (numeric)
#'     and `theta` (list), evaluating expected cumulative incidence
#'     at `times` (days since `date[1]`, as provided to [egf_init()])
#'     using parameter values `theta`. The default value of theta is
#'     `theta_mle`, hence the fitted cumulative curve is obtained as
#'     `eval_cum_inc(times)`.
#'   }
#'   \item{`eval_int_inc`}{A closure with arguments `time` (numeric)
#'     and `theta` (list), evaluating expected interval incidence.
#'     Returns `diff(eval_cum_inc(time, theta))`, hence the length
#'     of the value is `length(time)-1`. As with `eval_cum_inc`,
#'     the default value of theta is `theta_mle`.
#'   }
#'   \item{`init`}{Matches argument.}
#'   \item{`madf_out`}{The list output of [TMB::MakeADFun()].}
#'   \item{`optim_out`}{The list output of [stats::optim()].}
#'   \item{`call`}{The call to `egf()`, making the output
#'     reproducible with `eval(call)`.
#'   }
#' }
#'
#' @details
#' ## 1. Phenomenological models
#' Let \mjseqn{x(t)} be the expected number of cases observed
#' up to time \mjseqn{t} (i.e., expected cumulative incidence),
#' and let \mjseqn{x_0 = x(0) > 0}. Ignoring any baseline growth
#' (see Details 2), \mjseqn{x(t)} is modeled as an exponential,
#' logistic, or Richards curve.
#'
#' ### Exponential model
#' If \mjseqn{x(t)} follows
#'
#' \mjsdeqn{x'(t) = r x(t)\,,\qquad r > 0\,,}
#'
#' then \mjseqn{x(t)} grows exponentially as
#'
#' \mjsdeqn{x(t) = x_0 e^{r t}\,.}
#'
#' Hence the exponential model for \mjseqn{x(t)} requires
#' fitting two parameters to observed data: the exponential
#' growth rate \mjseqn{r} and initial value \mjseqn{x_0}.
#'
#' The exponential model ignores depletion of susceptible
#' individuals and implies continuous exponential growth
#' in \mjseqn{x(t)}. Hence it will only agree with epidemic
#' data during the (typically short) initial exponential
#' growth phase. \insertCite{Ma+14;textual}{epigrowthfit}
#' show that estimates of \mjseqn{r} are sensitive to the
#' choice of fitting window, and that more robust fits to
#' the data are likely to obtained with the logistic and
#' Richards models, at negligible cost.
#'
#' ### Logistic model
#' If \mjseqn{x(t)} follows
#'
#' \mjsdeqn{x'(t) = r x(t) \bigg(1 - \frac{x(t)}{K}\bigg)\,,\qquad r, K > 0\,,}
#'
#' and if \mjseqn{x_0 \in (0,K)}, then \mjseqn{x(t)} grows as
#'
#' \mjsdeqn{x(t) = \frac{K}{1 + \big(\frac{K}{x_0} - 1\big) e^{-r t}}}
#'
#' and increases to \mjseqn{K} as \mjseqn{t \to \infty}.
#' The logistic model can be reparametrized as
#'
#' \mjsdeqn{x(t) = \frac{K}{1 + e^{-r (t - t_\text{half})}}\,,}
#'
#' where \mjseqn{t_\text{half}} is the time at which
#' cumulative incidence attains half its final size,
#' satisfying \mjseqn{x(t_\text{half}) = \frac{K}{2}}.
#'
#' The reparametrized logistic model requires fitting
#' \mjseqn{r}, \mjseqn{K}, and \mjseqn{t_\text{half}}
#' to observed data.
#'
#' ### Richards
#' If \mjsdeqn{x(t)} follows
#'
#' \mjsdeqn{x'(t) = r x(t) \bigg(1 - \bigg(\frac{x(t)}{K}\bigg)^p\bigg)\,,\qquad r, K, p > 0\,,}
#'
#' and if \mjseqn{x_0 \in (0,K)}, then \mjseqn{x(t)} grows as
#'
#' \mjsdeqn{x(t) = \frac{K}{\big\lbrack 1 + \big(\big(\frac{K}{x_0}\big)^p - 1\big) e^{-r p t} \big\rbrack^{1/p}}}
#'
#' and increases to \mjseqn{K} as \mjseqn{t \to \infty}.
#' The Richards model can be reparametrized as
#'
#' \mjsdeqn{x(t) = \frac{K}{\big\lbrack 1 + (2^p - 1) * e^{-r p (t - t_\text{half}}) \big\rbrack^{1/p}}\,,}
#'
#' where, as with the logistic model, \mjseqn{t_\text{half}}
#' satisfies \mjseqn{x(t_\text{half}) = \frac{K}{2}}.
#'
#' The reparametrized logistic model requires fitting
#' \mjseqn{r}, \mjseqn{K}, \mjseqn{t_\text{half}}, and \mjseqn{p}
#' to observed data.
#'
#' ## 2. Baseline growth
#' For many historical epidemics, the observed data are counts
#' of all causes deaths, not counts of deaths due to the disease
#' of interest. Growth in disease mortality over time can still
#' be understood phenomenologically, provided that baseline
#' mortality (deaths unrelated to the epidemic) and disease
#' mortality are modeled separately. To account for baseline
#' mortality, it is assumed that deaths occur at a constant rate
#' \mjseqn{b > 0} in the absence of an epidemic. Then,
#' for example, the logistic model (see Details 1) becomes
#'
#' \mjsdeqn{x(t) = b t + \frac{K}{1 + \big(\frac{K}{x_0} - 1\big) e^{-r t}}\,,}
#'
#' where \mjseqn{x(t)} is to be interpreted as expected cumulative
#' mortality instead of expected cumulative incidence. Hence accounting
#' for baseline mortality requires that \mjseqn{b} is fit in addition
#' to the other model parameters.
#'
#' ## 3. Observation model
#' Let \mjseqn{Y(t_1,t_2)} be the number of cases observed between times
#' \mjseqn{t_1} and \mjseqn{t_2 > t_1} (i.e., observed interval incidence).
#' \mjseqn{Y(t_1,t_2)} is modeled as either a Poisson-distributed random
#' variable with mean \mjseqn{x(t_2) - x(t_1)},
#'
#' \mjsdeqn{Y(t_1,t_2) \sim \mathrm{Poisson}\big(x(t_2) - x(t_1)\big)\,,}
#'
#' or as negative binomial-distributed random variable with mean
#' \mjseqn{x(t_2) - x(t_1)} and dispersion \mjseqn{k > 0},
#'
#' \mjsdeqn{Y(t_1,t_2) \sim \mathrm{NegativeBinomial}\big(x(t_2) - x(t_1), k\big)\,.}
#'
#' The negative binomial observation model requires that \mjseqn{k} is fit
#' in addition to the other model parameters.
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
#' plot(x, inc = "interval")
#' plot(x, inc = "cumulative")
#'
#' @references
#' \insertRef{Ma+14}{epigrowthfit}
#'
#' @export
#' @importFrom TMB MakeADFun
#' @importFrom stats setNames optim
#' @useDynLib epigrowthfit
egf <- function(init, method = "BFGS", ...) {
  if (missing(init)) {
    stop("Missing argument `init`.")
  } else if (!inherits(init, "egf_init")) {
    stop("`init` must be an \"egf_init\" object.")
  }
  if (!is.character(method) || length(method) != 1 ||
        !method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-S", "SANN")) {
    stop("`method` must be an element of\n",
         "`c(\"Nelder-Mead\", \"BFGS\", \"CG\", \"L-BFGS-S\", \"SANN\")`.")
  }

  f <- init$first
  l <- init$last
  madf_data <- list(
    t = init$time[f:(l+1)],
    x = init$cases[f:l],
    curve_flag = match(init$curve, c("exponential", "logistic", "richards")) - 1,
    baseline_flag = 1 * init$include_baseline,
    distr_flag = match(init$distr, c("pois", "nbinom")) - 1
  )

  theta0 <- init$theta0
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

  theta_mle <- as.list(exp(optim_out$par))
  names(theta_mle) <- sub("log_", "", names(theta_mle))
  eval_cum_inc <- function(time, theta = theta_mle) {
    with(theta, {
      cum_inc <- switch(init$curve,
        exponential = x0 * exp(r * time),
        logistic = K / (1 + exp(-r * (time - thalf))),
        richards = K / (1 + (2^p - 1) * exp(-r * p * (time - thalf)))^(1 / p)
      )
      if (init$include_baseline) {
        cum_inc <- b * time + cum_inc
      }
      cum_inc
    })
  }
  eval_int_inc <- function(time, theta = theta_mle) {
    diff(eval_cum_inc(time, theta))
  }

  out <- list(
    theta_mle = as.list(theta_mle),
    eval_cum_inc = eval_cum_inc,
    eval_int_inc = eval_int_inc,
    init = init,
    madf_out = madf_out,
    optim_out = optim_out,
    call = match.call()
  )
  structure(out, class = c("egf", "list"))
}
