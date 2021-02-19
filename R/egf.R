#' Fit models of epidemic growth
#'
#' @description
#' Fits nonlinear mixed effects models of epidemic growth to one
#' or more disease incidence time series.
#'
#' @param formula_ts
#'   A formula of the form `y ~ x` or `y ~ x | group` specifying one
#'   or more incidence time series in long format, with `x` the name
#'   of a Date vector, `y` the name of a numeric vector, and `group`
#'   an interaction of one or more factors, splitting the data by
#'   time series.
#' @param formula_glmm
#'   A list of formulae of the form `par ~ terms`, where `par`
#'   is the name of a nonlinear model parameter (an element
#'   of `get_par_names(curve, distr, excess, link = TRUE)`)
#'   and `rhs` is a sum of mixed effects model terms following
#'   [`lme4`][lme4::lmer()]-like syntax. Alternatively,
#'   a formula of the form `~terms` to be recycled for all
#'   nonlinear model parameters.
#' @param data
#'   A data frame, list, or environment containing the variables
#'   named in `formula_ts` and `formula_glmm`.
#' @param window
#'   A factor of length `nrow(data)` such that `split(data, window)`
#'   splits `data` by fitting window, with `is.na(window)` indexing
#'   rows of `data` not belonging to a fitting window.
#' @param curve
#'   A character string specifying a cumulative incidence model.
#' @param distr
#'   A character string specifying an observation model.
#' @param excess
#'   A logical scalar. If `TRUE`, then a constant baseline
#'   mortality rate is estimated. Set to `TRUE` if what is
#'   observed (`y` if `formula = y ~ x`) is multiple causes
#'   mortality rather than disease mortality or disease incidence.
#' @param weekday
#'   A logical scalar. If `TRUE`, then day-of-week effects
#'   are estimated.
#' @param method
#'   A character string specifying an optimizer available through
#'   [stats::nlminb()], [stats::nlm()], or [stats::optim()].
#' @param na_action
#'   A character string indicating how `NA` in the incidence vector
#'   (`y` if `formula = y ~ x | group`) are handled. Note that `NA`
#'   in other variables are not tolerated.
#' @param sparse_X
#'   A logical scalar. If `TRUE`, then the fixed effects design matrix
#'   will be constructed in sparse format.
#' @param debug
#'   A logical scalar used for debugging. If `TRUE`, then `egf()`
#'   returns early with a list of optimization inputs.
#' @param par_init
#'   A full parameter vector for the first likelihood evaluation.
#'   Set to `NULL` to accept the internally generated default. If
#'   the default causes errors, then use `debug = TRUE` to retrieve
#'   its value, which will have the correct structure, and modify
#'   as necessary.
#' @param ...
#'   Optional arguments to the optimizer specified by `method`.
#'
#' @return
#' If `debug = FALSE`, then a list inheriting from class `"egf"`,
#' with elements:
#' \item{`frame`}{
#'   The model frame. See [make_frame()].
#' }
#' \item{`window`}{
#'   A factor of length `nrow(data)` such that `split(frame, window)`
#'   splits the model frame by fitting window.
#' }
#' \item{`curve`, `distr`, `excess`, `weekday`, `method`}{
#'   Copies of the so-named arguments (after matching).
#' }
#' \item{`tmb_args`}{
#'   A list of arguments to [TMB::MakeADFun()]. See [make_tmb_args()].
#' }
#' \item{`tmb_out`}{
#'   The list output of [TMB::MakeADFun()] (after optimization).
#' }
#' \item{`optim_out`}{
#'   The list output of the optimizer specified by `method`.
#' }
#' \item{`par_init`}{
#'   The full parameter vector of the first likelihood evaluation.
#'   Matches the so-named argument if given.
#' }
#' \item{`par`}{
#'   The full parameter vector of the best likelihood evaluation.
#' }
#' \item{`nonrandom`}{
#'   An integer vector indexing the nonrandom segment of `par`.
#' }
#' \item{`report`}{
#'   A list containing all variables `REPORT()`ed and `ADREPORT()`ed
#'   in the C++ template, with an additional element `cov` giving the
#'   covariance matrix corresponding to `par[nonrandom]`.
#' }
#' \item{`nll`}{
#'   A numeric scalar giving the value of the negative log
#'   Laplace approximation of the marginal likelihood function
#'   at `par[nonrandom]`.
#' }
#' \item{`nll_func`, `nll_grad`}{
#'   Closures taking a numeric vector `x` as an argument and returning
#'   the value or gradient of the negative log Laplace approximation of
#'   the marginal likelihood function at `x[nonrandom]`.
#' }
#' \item{`call`}{
#'   The call to `egf()`, allowing for updates to the `"egf"` object
#'   via [stats::update()].
#' }
#' If `debug = TRUE`, then a list with only the elements `frame`,
#' `tmb_args`, and `par_init`.
#'
#' @export
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats nlminb nlm optim
#' @useDynLib epigrowthfit
egf <- function(formula_ts,
                formula_glmm,
                data = parent.frame(),
                window,
                curve = c("logistic", "richards", "exponential", "subexponential", "gompertz"),
                distr = c("nbinom", "pois"),
                excess = FALSE,
                weekday = FALSE,
                method = c("nlminb", "nlm", "Nelder-Mead", "BFGS"),
                na_action = c("pass", "fail"),
                sparse_X = FALSE,
                debug = FALSE,
                par_init = NULL,
                ...) {
  curve <- match.arg(curve)
  distr <- match.arg(distr)
  stop_if_not_tf(excess)
  stop_if_not_tf(weekday)
  method <- match.arg(method)
  na_action <- match.arg(na_action)
  stop_if_not_tf(sparse_X)
  stop_if_not_tf(debug)

  frame <- make_frame(
    formula_ts = formula_ts,
    formula_glmm = formula_glmm,
    data = data,
    window = window,
    curve = curve,
    distr = distr,
    excess = excess,
    weekday = weekday,
    na_action = na_action
  )
  tmb_args <- make_tmb_args(
    frame = frame,
    curve = curve,
    distr = distr,
    excess = excess,
    sparse_X = sparse_X,
    par_init = par_init
  )

  if (debug) {
    pl <- tmb_args$parameters
    for (s in names(pl)) {
      names(pl[[s]]) <- rep.int(s, length(pl[[s]]))
    }
    out <- list(
      frame = frame,
      tmb_args = tmb_args,
      par_init = rename_par(unlist(pl))
    )
    return(out)
  }

  tmb_out <- do.call(MakeADFun, tmb_args)
  optim_out <- optim_tmb_out(tmb_out, method = method, ...)

  par_init <- rename_par(tmb_out$env$par)
  par <- rename_par(tmb_out$env$last.par.best)
  nonrandom <- grep("^b\\[", names(par), invert = TRUE)

  s <- switch(method, nlminb = "objective", nlm = "minimum", "value")
  nll <- optim_out[[s]]
  nll_func <- function(x = par) as.numeric(tmb_out$fn(x[nonrandom]))
  nll_grad <- function(x = par) as.vector(tmb_out$gr(x[nonrandom]))

  sdr <- sdreport(tmb_out)
  report <- c(
    tmb_out$report(par),
    list(cov = sdr$cov.fixed),
    split_sdreport(sdr)
  )

  out <- list(
    frame = frame,
    index = attr(frame, "index"),
    curve = curve,
    distr = distr,
    excess = excess,
    method = method,
    tmb_args = tmb_args,
    tmb_out = tmb_out,
    optim_out = optim_out,
    par_init = par_init,
    par = par,
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
