#' Fit models of epidemic growth
#'
#' @description
#' Fits nonlinear mixed effects models of epidemic growth to one
#' or more disease incidence time series.
#'
#' @param formula_ts
#'   A formula of the form `x ~ time` or `x ~ time | ts` specifying one
#'   or more incidence time series in long format. `time` must evaluate
#'   to a numeric or Date vector. Numeric `time` is assumed to measure
#'   time as a number of days since _the end of_ date `origin`. Date
#'   `time` is coerced to numeric `julian(time, origin)` (see Details).
#'   `x` must evaluate to a numeric vector. Within a time series,
#'   `x[i]` should specify the number of cases observed from `time[i-1]`
#'   to `time[i]`. Finally, `ts` must evaluate to a factor, such that
#'   `split(data.frame(time, x), ts)` returns a list of time series.
#'   Note that `x ~ time` is equivalent to `x ~ time | ts` with `ts`
#'   set equal to `rep(factor(1), length(x))`.
#' @param formula_par
#'   A named list of formulae of the form `~terms`, specifying
#'   mixed effects models ([`lme4`][lme4::lmer()]-like syntax)
#'   for nonlinear model parameters. `names(formula_par)` must
#'    be a subset of
#'   `get_par_names(curve, distr, excess, weekday, link = TRUE)`.
#'   `~1` is the default for parameters not assigned a formula.
#'   Alternatively, `formula_par` may itself be a formula of
#'   the form `~terms`. In this case, the formula is recycled
#'   for all nonlinear model parameters. Mixed effects variables
#'   must match the length of the time series variables and be
#'   constant in each level of `window`.
#' @param data
#'   A data frame, list, or environment. `egf()` looks here
#'   for variables used in `formula_ts` and `formula_par`,
#'   before looking in formula environments.
#' @param window
#'   A factor of length `nrow(data)` such that `split(data, window)`
#'   splits `data` by fitting window, with `is.na(window)` indexing
#'   rows of `data` not belonging to a fitting window. Note that
#'   such a factor can be constructed from a table of fitting window
#'   endpoints using [make_window()].)
#' @param origin
#'   A Date specifying a reference time (see `formula_ts`).
#' @param curve
#'   A character string specifying a cumulative incidence model.
#' @param excess
#'   A logical scalar. If `TRUE`, then a constant baseline mortality
#'   rate is estimated. Set to `TRUE` if what is observed is
#'   multiple causes mortality rather than disease mortality
#'   or disease incidence.
#' @param distr
#'   A character string specifying an observation model.
#' @param weekday
#'   An integer or logical scalar. If `weekday > 0`, then weekday
#'   effects are estimated as offsets relative to the indicated day
#'   (Sunday if `weekday = 1`, Monday if `weekday = 2`, and so on).
#'   Currently, weekday effect estimation requires time (that is,
#'   the value of `foo` in `formula_ts = x ~ foo | ts`) to be an
#'   integer (in the sense of `all.equal(foo, round(foo))`) or Date
#'   vector with 1-day spacing in all fitting windows.
#' @param method
#'   A character string specifying an optimizer available through
#'   [stats::nlminb()], [stats::nlm()], or [stats::optim()].
#' @param na_action
#'   A character string indicating how `NA` in incidence (`x`
#'   if `formula_ts = x ~ time | ts`) are handled. Note that `NA`
#'   in other `formula_ts` variables are an error, and `NA` in
#'   `formula_par` variables within fitting windows are an error.
#' @param sparse_X
#'   A logical scalar. If `TRUE`, then the fixed effects design
#'   matrix is constructed in sparse format.
#' @param debug
#'   A logical scalar used for debugging. If `TRUE`, then
#'   `egf()` returns early with a list of optimization inputs.
#' @param init
#'   A full parameter vector for the first likelihood evaluation.
#'   Set to `NULL` to accept the internally generated default.
#'   Use `debug = TRUE` to retrieve this default, which is useful
#'   as a template.
#' @param ...
#'   Optional arguments to the optimizer specified by `method`.
#'
#' @details
#' Coercion of Date `time` to numeric `julian(time, origin)` assumes
#' that each element `time[i]` can be read as "end of date `time[i]`",
#' so that observation `x[i]` in a given time series is the number of
#' cases observed from the end of date `time[i-1]` to the end of date
#' `time[i]`.
#'
#' @return
#' If `debug = FALSE`, then a list inheriting from class `"egf"`,
#' with elements:
#' \item{`frame_ts`}{
#'   The time series model frame, constructed from `formula_ts`.
#'   `origin` is retained as an attribute. See [make_frames()].
#' }
#' \item{`frame_par`}{
#'   A list of mixed effects model frames (one per nonlinear model
#'   parameter), constructed from `formula_par` (after completion).
#'   See [make_frames()].
#' }
#' \item{`endpoints`}{
#'   A data frame with variables `ts`, `window`, `start`, and `end`
#'   listing fitting window endpoints as numbers of days since time
#'   23:59:59 on a reference date, which is stored as an attribute.
#' }
#' \item{`curve`, `excess`, `distr`, `weekday`, `method`}{
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
#' \item{`init`}{
#'   The full parameter vector of the first likelihood evaluation.
#' }
#' \item{`best`}{
#'   The full parameter vector of the best likelihood evaluation.
#' }
#' \item{`nonrandom`}{
#'   An integer vector indexing the nonrandom segment of `best`.
#' }
#' \item{`report`}{
#'   A list containing all variables `REPORT()`ed and `ADREPORT()`ed
#'   in the C++ template, with an additional element `cov` giving the
#'   covariance matrix corresponding to `par[nonrandom]`.
#' }
#' \item{`nll`}{
#'   A numeric scalar giving the value of the negative log
#'   Laplace approximation of the marginal likelihood function
#'   at `best[nonrandom]`.
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
#' If `debug = TRUE`, then a list with only the elements `frame_ts`,
#' `frame_par`, `init`, and `tmb_args`.
#'
#' @export
#' @importFrom TMB MakeADFun sdreport
#' @useDynLib epigrowthfit
egf <- function(formula_ts,
                formula_par,
                data = parent.frame(),
                window,
                origin = .Date(0),
                curve = c("logistic", "richards", "exponential", "subexponential", "gompertz"),
                excess = FALSE,
                distr = c("nbinom", "pois"),
                weekday = FALSE,
                method = c("nlminb", "nlm", "Nelder-Mead", "BFGS"),
                na_action = c("pass", "fail"),
                sparse_X = FALSE,
                debug = FALSE,
                init = NULL,
                ...) {
  stop_if_not(
    inherits(origin, "Date"),
    length(origin) == 1L,
    !is.na(origin),
    m = "`origin` must be a Date vector of length 1."
  )
  curve <- match.arg(curve)
  distr <- match.arg(distr)
  stop_if_not_true_false(excess)
  if (is.logical(weekday) && length(weekday) == 1L && !is.na(weekday)) {
    weekday <- as.integer(weekday)
  } else {
    stop_if_not_integer(weekday)
    ## Coercion to element of 0:7
    weekday <- (weekday > 0) * as.integer(1 + (weekday - 1) %% 7)
  }
  method <- match.arg(method)
  na_action <- match.arg(na_action)
  stop_if_not_true_false(sparse_X)
  stop_if_not_true_false(debug)

  frames <- make_frames(
    formula_ts = formula_ts,
    formula_par = formula_par,
    data = data,
    window = window,
    origin = origin,
    curve = curve,
    distr = distr,
    excess = excess,
    weekday = weekday,
    na_action = na_action,
    init = init
  )
  tmb_args <- make_tmb_args(
    frame_ts = frames$frame_ts,
    frame_par = frames$frame_par,
    curve = curve,
    distr = distr,
    excess = excess,
    weekday = weekday,
    sparse_X = sparse_X,
    init = init
  )

  if (debug) {
    init_split <- tmb_args$parameters
    if (!has_random(tmb_args$data)) {
      init_split <- init_split["beta"]
    }
    init <- unlist(init_split)
    names(init) <- enum_dupl_string(rep.int(names(init_split), lengths(init_split)))
    out <- list(
      frame_ts = frames$frame_ts,
      frame_par = frames$frame_par,
      init = init,
      tmb_args = tmb_args
    )
    return(out)
  }

  tmb_out <- do.call(MakeADFun, tmb_args)
  optim_out <- optim_tmb_out(tmb_out, method = method, ...)

  init <- enum_dupl_names(tmb_out$env$par)
  best <- enum_dupl_names(tmb_out$env$last.par.best)
  nonrandom <- grep("^b\\[", names(best), invert = TRUE)

  s <- switch(method, nlminb = "objective", nlm = "minimum", "value")
  nll <- optim_out[[s]]
  nll_func <- function(x = par) as.numeric(tmb_out$fn(x[nonrandom]))
  nll_grad <- function(x = par) as.numeric(tmb_out$gr(x[nonrandom]))

  sdr <- sdreport(tmb_out)
  report <- c(
    tmb_out$report(best),
    list(cov = sdr$cov.fixed),
    split_sdreport(sdr)
  )

  out <- list(
    frame_ts = frames$frame_ts,
    frame_par = frames$frame_par,
    endpoints = get_endpoints(frames$frame_ts),
    curve = curve,
    excess = excess,
    distr = distr,
    weekday = weekday,
    method = method,
    tmb_args = tmb_args,
    tmb_out = tmb_out,
    optim_out = optim_out,
    init = init,
    best = best,
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
