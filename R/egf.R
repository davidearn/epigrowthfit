#' \loadmathjax
#' Fit a model of epidemic growth
#'
#' @description
#' Fits phenomenological models of epidemic growth to one or more
#' disease incidence time series. Accounts for variation of model
#' parameters (e.g., the initial exponential growth rate \mjseqn{r})
#' across time series (e.g., between epidemic waves or geographical
#' units) by assuming that the parameters follow a user-specified
#' generalized linear mixed effects model. First and second order
#' derivatives of the likelihood are calculated
#' by automatic differentiation using package **TMB**.
#'
#' @param formula
#'   A formula of the form `y ~ x` specifying an incidence time series
#'   or multiple incidence time series in long format. `x` and `y` must
#'   be vectors of dates (see `date_format`) and integers, respectively.
#' @param fixed,random
#'   Named lists of formulae of the form `~rhs`, together specifying
#'   a generalized linear mixed effects model (the fixed and random
#'   components) for each parameter of the _incidence_ model indicated
#'   by `curve`, `distr`, and `excess`. A list of valid names can be
#'   obtained with `get_par_names(curve, distr, excess, link = TRUE)`.
#'   Alternatively, formulae (rather than lists of formulae) to be
#'   recycled for all parameters. Syntax is [`lme4`][lme4::lmer()]-like,
#'   but here all variables must be factors. `fixed` formulae are
#'   restricted to one term and so must have the form `~f1:...:fn`
#'   or `~1` (the default for each missing list element). `random`
#'   formulae must give sums of terms of the form `(1 | f1:...:fn)`,
#'   `(1 | f1/.../fn)`, or `(1 | f1 * ... * fn)`. Use `NULL` (the
#'   default for each missing list element) instead of a formula to
#'   indicate absence of random effects.
#' @param data
#'   A data frame, list, or environment containing the variables
#'   named in `formula`, `fixed`, and `random`.
#' @param index
#'   A factor of length `nrow(data)` such that `split(data, index)`
#'   splits `data` by fitting window. In this case, `is.na(index)`
#'   should index rows of `data` not belonging to a fitting window.
#'   `NULL` (the default) is equivalent to `factor(integer(nrow(data)))`.
#' @param curve
#'   A character string specifying a cumulative incidence model.
#' @param distr
#'   A character string specifying an observation model.
#' @param excess
#'   A logical scalar. If `TRUE`, then a constant baseline mortality
#'   rate will be estimated. Set to `TRUE` if what is observed
#'   (`y` if `formula = y ~ x`) is multiple causes mortality
#'   rather than disease mortality or incidence.
#' @param method
#'   A character string specifying an optimizer available through
#'   [stats::nlminb()], [stats::nlm()], or [stats::optim()].
#' @param na_action
#'   A character string indicating how `NA` in the incidence vector
#'   (`y` if `formula = y ~ x`) are handled. See [stats::na.fail()].
#' @param sparse_X
#'   A logical scalar. If `TRUE`, then the fixed effects design matrix
#'   will be constructed in sparse format.
#' @param date_format
#'   A character vector passed to [as.Date()] argument `tryFormats`,
#'   used to convert dates (`x` if `formula = y ~ x`) from character
#'   to Date, if necessary. See [base::strptime()] for a list of
#'   conversion specifications.
#' @param debug
#'   A logical scalar used for debugging. If `TRUE`, then `egf()`
#'   returns early with a list of optimization inputs.
#' @param par_init
#'   A full parameter vector for the first likelihood evaluation.
#'   Set to `NULL` to accept the internally generated default. If
#'   the default causes errors, then use `debug = TRUE` to retrieve
#'   its value, which will have the correct structure, and modify as
#'   desired.
#' @param ...
#'   Optional arguments to the optimizer specified by `method`.
#'
#' @return
#' If `debug = FALSE`, then an `"egf"` object, which is a list with
#' elements:
#' \item{`frame`}{
#'   The model frame. This is a data frame containing _only_ the
#'   variables named in `formula`, `fixed`, and `random`. Data not
#'   belonging to a fitting window appear in the tail.
#' }
#' \item{`index`}{
#'   A factor of length `nrow(data)` such that `split(frame, index)`
#'   splits the model frame by fitting window. Not necessarily identical
#'   to the so-named argument, as `frame` may be a permutation of the
#'   rows of `data`.
#' }
#' \item{`curve`, `distr`, `excess`}{
#'   Copies of the so-named arguments.
#' }
#' \item{`tmb_args`}{
#'   A list of arguments to [TMB::MakeADFun()]. See [make_tmb_args()].
#' }
#' \item{`tmb_out`}{
#'   The (optimized) list output of [TMB::MakeADFun()].
#' }
#' \item{`optim_out`}{
#'   The list output of the optimizer.
#' }
#' \item{`par_init`}{
#'   The full parameter vector of the first likelihood evaluation.
#'   Matches the so-named argument if given.
#' }
#' \item{`par`}{
#'   The full parameter vector of the best likelihood evaluation.
#' }
#' \item{`par_info`}{
#'   A data frame summarizing the meaning of the elements of `par`
#'   in relation to the mixed effects model. See [get_par_info()].
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
#'   A numeric scalar giving the negative log Laplace approximation
#'   of the marginal likelihood of `par[nonrandom]`.
#' }
#' \item{`nll_func`, `nll_grad`}{
#'   Closures taking a numeric vector `x` of length `length(nonrandom)`
#'   as an argument, and returning the value of the negative log
#'   Laplace approximation of the marginal likelihood and its gradient,
#'   respectively, evaluated at `x`.
#' }
#' \item{`call`}{
#'   The call to `egf()`, allowing for updates to the `"egf"` object
#'   via [stats::update()].
#' }
#'
#' If `debug = TRUE`, then a list containing `frame`, `tmb_args`,
#' `par_init`, and `par_init_info`.
#'
#' @export
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats nlminb nlm optim
egf <- function(formula,
                fixed = ~1,
                random = NULL,
                data = parent.frame(),
                index = NULL,
                curve = c("logistic", "richards", "exponential", "subexponential", "gompertz"),
                distr = c("nbinom", "pois"),
                excess = FALSE,
                method = c("nlminb", "nlm", "Nelder-Mead", "BFGS"),
                na_action = c("pass", "fail"),
                sparse_X = FALSE,
                date_format = "%Y-%m-%d",
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

  frame <- make_frame(formula, fixed, random, data, index,
                      curve, distr, excess, na_action, date_format)
  tmb_args <- make_tmb_args(frame, curve, distr, excess,
                            sparse_X, par_init)
  if (debug) {
    pl <- tmb_args$parameters
    par_init <- rename_par(unlist(Map("names<-", unname(pl), Map(rep, names(pl), lengths(pl)))))
    par_init_info <- get_par_info(par_init, tmb_args$data)
    out <- list(
      frame = frame,
      tmb_args = tmb_args,
      par_init = par_init,
      par_init_info = par_init_info
    )
    return(out)
  }

  tmb_out <- do.call(MakeADFun, tmb_args)
  optim_out <- optim_tmb_out(tmb_out, method, ...)
  nll_name <- switch(method, nlminb = "objective", nlm = "minimum", "value")
  nll <- optim_out[[nll_name]]

  par_init <- rename_par(tmb_out$env$par)
  par <- rename_par(tmb_out$env$last.par.best)
  par_info <- get_par_info(par, tmb_args$data)
  nonrandom <- grep("^b\\[", names(par), invert = TRUE)

  ## Store only necessary variables in the function environment
  ## to minimize duplication
  nll_func <- evalq(
    expr = function(x = par) as.numeric(fn(x[nonrandom])),
    envir = list(par = par, nonrandom = nonrandom, fn = tmb_out$fn),
    enclos = baseenv()
  )
  nll_grad <- evalq(
    expr = function(x = par) as.vector(gr(x[nonrandom])),
    envir = list(par = par, nonrandom = nonrandom, gr = tmb_out$gr),
    enclos = baseenv()
  )

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
    tmb_args = tmb_args,
    tmb_out = tmb_out,
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
