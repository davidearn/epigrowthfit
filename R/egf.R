#' Fit models of epidemic growth
#'
#' @description
#' Fits nonlinear mixed effects models of epidemic growth to one
#' or more disease incidence time series. First and second order
#' derivatives of the negative log likelihood (more precisely,
#' the negative log Laplace approximation of the marginal likelihood)
#' are calculated by automatic differentiation using R package **TMB**.
#'
#' @param formula
#'   A formula of the form `y ~ x` specifying one or more incidence
#'   time series in long format, with `x` Date and `y` numeric.
#' @param fixed,random
#'   Named lists of formulae of the form `~rhs`, together specifying
#'   a generalized linear mixed effects model (the fixed and random
#'   components) for each parameter of the _incidence_ model indicated
#'   by `curve`, `distr`, and `excess`.
#'   A list of valid names can be obtained with
#'   `get_par_names(curve, distr, excess, link = TRUE)`.
#'   Alternatively, formulae (rather than lists of formulae) to be
#'   recycled for all parameters.
#'   Syntax is [`lme4`][lme4::lmer()]-like, but here all variables
#'   must be factors.
#'   `fixed` formulae are restricted to one term and so must be `~1`
#'   (the default for each missing list element) or have the form
#'   `~f1:...:fn`.
#'   `random` formulae are restricted to sums of terms of the form
#'   `(1 | f1:...:fn)`, `(1 | f1/.../fn)`, or `(1 | f1 * ... * fn)`.
#'   Use `NULL` (the default for each missing list element) instead
#'   of a formula to indicate absence of random effects.
#' @param group_by
#'   A formula of the form `~f1:...:fn`,
#'   such that `split(data, interaction(data[all.vars(group_by)]))`
#'   splits `data` by time series. Use `~1` (the default) if there
#'   is only one time series.
#' @param data
#'   A data frame, list, or environment containing the variables
#'   named in `formula`, `fixed`, `random`, and `group_by`. Missing
#'   values in incidence (`y` if `formula = y ~ x`) are tolerated
#'   only if `na_action = "pass"`. Missing values in other variables
#'   are not tolerated.
#' @param index
#'   A factor of length `nrow(data)` such that `split(data, index)`
#'   splits `data` by fitting window, with `is.na(index)` indexing
#'   rows of `data` not belonging to a fitting window. `NULL`
#'   (the default) is equivalent to `rep(factor(0), nrow(data))`.
#' @param curve
#'   A character string specifying a cumulative incidence model.
#' @param distr
#'   A character string specifying an observation model.
#' @param excess
#'   A logical scalar. If `TRUE`, then a constant baseline
#'   mortality rate is estimated. Set to `TRUE` if what is
#'   observed (`y` if `formula = y ~ x`) is multiple causes
#'   mortality rather than disease mortality or disease incidence.
#' @param method
#'   A character string specifying an optimizer available through
#'   [stats::nlminb()], [stats::nlm()], or [stats::optim()].
#' @param na_action
#'   A character string indicating how `NA` in the incidence vector
#'   (`y` if `formula = y ~ x`) are handled. See [stats::na.fail()].
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
#' \item{`index`}{
#'   A factor of length `nrow(data)` such that `split(frame, index)`
#'   splits the model frame by fitting window. Not usually identical
#'   to the so-named argument.
#' }
#' \item{`curve`, `distr`, `excess`, `method`}{
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
egf <- function(formula,
                fixed = ~1,
                random = NULL,
                group_by = ~1,
                data = parent.frame(),
                index = NULL,
                curve = c("logistic", "richards", "exponential", "subexponential", "gompertz"),
                distr = c("nbinom", "pois"),
                excess = FALSE,
                method = c("nlminb", "nlm", "Nelder-Mead", "BFGS"),
                na_action = c("pass", "fail"),
                sparse_X = FALSE,
                debug = FALSE,
                par_init = NULL,
                ...) {
  curve <- match.arg(curve)
  distr <- match.arg(distr)
  stop_if_not_tf(excess)
  method <- match.arg(method)
  na_action <- match.arg(na_action)
  stop_if_not_tf(sparse_X)
  stop_if_not_tf(debug)

  frame <- make_frame(
    formula = formula,
    fixed = fixed,
    random = random,
    group_by = group_by,
    data = data,
    index = index,
    curve = curve,
    distr = distr,
    excess = excess,
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
