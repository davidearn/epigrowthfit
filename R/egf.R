#' \loadmathjax
#' Fit a model of epidemic growth
#'
#' @description
#' Fits phenomenological models of epidemic growth to one or more
#' disease incidence time series. Accounts for variation of model
#' parameters (e.g., the initial exponential growth rate \mjseqn{r})
#' across time series (e.g., between epidemic waves or geographical
#' units) by assuming a user-specified generalized linear mixed model
#' for the parameters. First and second order derivatives of the
#' likelihood are calculated by automatic differentiation using
#' package **TMB**.
#'
#' @param formula
#'   A formula of the form `y ~ x` locating incidence time series
#'   in `data`.
#' @param fixed,random
#'   A formula or named list of formula of form `~rhs` specifying
#'   fixed and random effects. See Details.
#' @param data
#'   A data frame, list, or environment containing the variables
#'   in `formula`, `fixed`, and `random`.
#' @param index
#'   A factor of length `nrow(data)` specifying fitting windows
#'   (subsets of data )
#'
#' @details
#' ## 1. Formula specification
#' For fixed effects, `rhs` must be a
#'   single interaction of `n >= 1` factors, as in `f1:...:fn`.
#'   For random effects, `rhs` must be a sum of terms of the form

#' @importFrom TMB MakeADFun sdreport
egf <- function(formula,
                fixed = ~1,
                random = NULL,
                data = parent.frame(),
                index = NULL,
                curve = c("logistic", "richards", "exponential", "subexponential", "gompertz"),
                distr = c("nbinom", "pois"),
                excess = FALSE,
                method = c("nlminb", "nlm", "Nelder-Mead", "BFGS", "CG"),
                na_action = c("exclude", "fail"),
                sparse_X = FALSE,
                date_format = "%Y-%m-%d",
                ...) {
  curve <- match.arg(curve)
  distr <- match.arg(distr)
  stop_if_not(
    inherits(excess, "logical"),
    length(excess) == 1L,
    !is.na(excess),
    m = "`excess` must be TRUE or FALSE."
  )
  na_action <- match.arg(na_action)
  method <- match.arg(method)

  formula <- check_formula(formula)
  fixed <- check_fixed(fixed, curve, distr, excess)
  random <- check_random(random, curve, distr, excess)
  frame <- check_data(formula, fixed, random, data, index,
                      na_action, date_format)

  madf_data <- make_madf_data(frame, curve, distr, excess, sparse_X)
  madf_args <- list(
    data = madf_data,
    parameters = make_madf_parameters(madf_data, curve),
    map = make_madf_map(madf_data),
    random = make_madf_random(madf_data),
    DLL = "epigrowthfit",
    silent = TRUE
  )
  madf_out <- do.call(MakeADFun, madf_args)

  if (method == "nlminb") {
    optim_out <- nlminb(
      start = madf_out$par,
      objective = madf_out$fn,
      gradient = madf_out$gr,
      ...
    )
    nll <- optim_out$objective
  } else if (method == "nlm") {
    optim_out <- nlm(
      f = structure(madf_out$fn, gradient = madf_out$gr),
      p = madf_out$par,
      ...
    )
    nll <- optim_out$minimum
  } else {
    optim_out <- optim(
      par = madf_out$par,
      fn = madf_out$fn,
      gr = madf_out$gr,
      method = method,
      ...
    )
    nll <- optim_out$value
  }

  par_init <- rename_par(madf_out$env$par)
  par <- rename_par(madf_out$env$last.par.best)
  par_info <- get_par_info(par = par, madf_data = madf_data)
  nonrandom <- grep("^b\\[", names(par), invert = TRUE)

  ## Store only necessary variables in environment(f)
  nll_func <- evalq(
    expr = function(x = par) as.numeric(fn(x[nonrandom])),
    envir = list(par = par, nonrandom = nonrandom, fn = madf_out$fn),
    enclos = baseenv()
  )
  nll_grad <- evalq(
    expr = function(x = par) as.vector(gr(x[nonrandom])),
    envir = list(par = par, nonrandom = nonrandom, gr = madf_out$gr),
    enclos = baseenv()
  )

  sdr <- sdreport(madf_out)
  report <- c(
    madf_out$report(par),
    list(cov = sdr$cov.fixed),
    split_sdreport(sdr)
  )

  out <- list(
    frame = frame,
    index = attr(frame, "index"),
    curve = curve,
    distr = distr,
    excess = excess,
    madf_args = madf_args,
    madf_out = madf_out,
    optim_out = optim_out,
    par_init = par_init,
    par = par,
    par_info = par_info,
    nonrandom = nonrandom,
    report = report,
    nll = nll,
    nll_func = nll_func,
    nll_grad = nll_grad,
    call = match.call()
  )
  class(out) <- c("egf", "list")
  out
}
