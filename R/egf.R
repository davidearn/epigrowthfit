#' Fit models of epidemic growth
#'
#' @description
#' Fits nonlinear mixed effects models of epidemic growth to one
#' or more disease incidence time series.
#'
#' @param formula_ts
#'   A formula of the form `x ~ time` or `x ~ time | ts` specifying
#'   one or more incidence time series in long format. `time` must
#'   evaluate to a numeric or Date vector. Numeric `time` is assumed
#'   to measure time as a number of days since _the end of_ date
#'   `origin`. Date `time` is coerced to numeric `julian(time, origin)`
#'   (see Details). `x` must evaluate to a numeric vector. Within a
#'   time series, `x[i]` should specify the number of cases observed
#'   from `time[i-1]` to `time[i]`. Finally, `ts` must evaluate to
#'   a factor such that `split(data.frame(time, x), ts)` returns
#'   a list of time series. Note that `x ~ time` is equivalent to
#'   `x ~ time | ts` with `ts` set equal to `rep(factor(1), length(x))`.
#' @param formula_par
#'   A list of formulae of the form `par ~ terms` specifying
#'   mixed effects models ([`lme4`][lme4::lmer()]-like syntax)
#'   for nonlinear model parameters. `deparse(par)` must be an
#'   element of
#'   `get_par_names(curve, distr, excess, weekday, link = TRUE)`.
#'   `~1` is the default for parameters not assigned a formula.
#'   Alternatively, `formula_par` may be a formula of the form
#'   `~terms`. In this case, the formula is recycled for all
#'   nonlinear model parameters. Note that "individuals" in each
#'   model are fitting windows, and model frames constructed
#'   from `formula_par` and `data_par` are expected to correspond
#'   rowwise to `endpoints` (see Details).
#' @param data_ts,data_par
#'   Data frames, lists, or environments. These are searched prior
#'   to formula environments for variables used in `formula_ts` and
#'   `formula_par`, respectively.
#' @param endpoints
#'   A data frame, list, or environment with variables `start`
#'   and `end`, and any further variables necessary to evaluate
#'   `ts` if `formula_ts = x ~ time | ts`. `start` and `end`
#'   must be numeric or Date vectors listing start and end times
#'   for all fitting windows. `ts` must evaluate to a factor
#'   indicating the time series in which each window is found.
#'   Within time series, intervals `[start[i], end[i]]` must be
#'   disjoint and contain at least two time points from `time`.
#' @param origin
#'   A Date specifying a reference time.
#' @param curve
#'   A character string specifying a cumulative incidence model.
#' @param excess
#'   A logical scalar. If `TRUE`, then a constant baseline mortality
#'   rate is estimated. Set to `TRUE` if what is observed is multiple
#'   causes mortality rather than disease mortality or disease incidence.
#' @param distr
#'   A character string specifying an observation model.
#' @param weekday
#'   An integer or logical scalar. If `weekday > 0`, then weekday
#'   effects are estimated as offsets relative to the indicated day
#'   (Sunday if `weekday = 1`, Monday if `weekday = 2`, and so on).
#'   Currently, weekday effect estimation requires `time` in
#'   `formula_ts = x ~ time | ts` to evaluate to an integer
#'   (in the sense of `all.equal(time, round(time))`) or Date vector
#'   with 1-day spacing in all fitting windows.
#' @param method
#'   A character string specifying an optimizer available through
#'   [stats::nlminb()], [stats::nlm()], or [stats::optim()].
#' @param na_action
#'   A character vector to be recycled to length 2.
#'   `na_action[1L]` affects the handling of `NA`
#'   in `x` if `formula_ts = x ~ time | ts`.
#'   `"fail"` is to throw an error.
#'   `"exclude"` is to ignore `NA` when fitting and retain `NA`
#'   in predictions.
#'   `"pass"` is to ignore `NA` when fitting and predict `NA`.
#'   `NA` in `time` and `ts` are an error regardless of `na_action`.
#'   `na_action[2L]` affects the handling of `NA`
#'   in `formula_par` variables.
#'   `"fail"` is to throw an error.
#'   `"exclude"` and `"pass"` are to discard fitting windows
#'   with incomplete data.
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
#' @param append
#'   An expression indicating variables in `data_par` to be preserved
#'   in the returned `"egf"` object for use by methods.
#'   It is evaluated similarly to the `select` argument of [subset()].
#'   Currently, usage is supported only if `data_par` is a data frame.
#' @param ...
#'   Optional arguments to the optimizer specified by `method`.
#'
#' @details
#' If `formula_ts = x ~ time | ts`, then coercion of Date `time` to
#' numeric `julian(time, origin)` assumes that each element `time[i]`
#' can be read as "end of date `time[i]`", so that observation `x[i]`
#' in a given time series is the number of cases observed from the
#' end of date `time[i-1]` to the end of date `time[i]`.
#'
#' To avoid unexpected mismatch between `endpoints` and mixed effects
#' model frames constructed from `formula_par` and `data_par`, it is
#' helpful to keep all `endpoints` and `formula_par` variables in a
#' common data frame `d` and set `endpoints = d` and `data_par = d`.
#'
#' @return
#' If `debug = FALSE`, then a list inheriting from class `"egf"`,
#' with elements:
#' \item{`endpoints`}{
#'   A data frame with variables `ts`, `window`, `start`, and
#'   `end` listing start and end times for all fitting windows.
#'   (Supplied intervals are contracted internally so that
#'   `start` and `end` are precisely the minimum and maximum
#'   of the set of time points contained in the window.)
#'   Rows are ordered by time series and chronologically
#'   within time series. `origin` is retained as an attribute.
#' }
#' \item{`frame_ts`}{
#'   The time series model frame, constructed from `formula_ts`
#'   and `data_ts`, with variables `ts`, `window`, `time`, and
#'   `x`. `ts` and `window` are factors that can be used to split
#'   `frame_ts` by time series and by fitting window, respectively.
#'   Rows are ordered by time series and chronologically within
#'   time series. `terms(formula_ts)` and `origin` are retained
#'   as attributes. `unclass(window)` indexes rows of `endpoints`.
#'   That is, time series data for the fitting window defined by
#'   row `i` of `endpoints` can be found in the rows of `frame_ts`
#'   for which `window = levels(window)[i]`.
#' }
#' \item{`frame_par`}{
#'   A list of mixed effects model frames, constructed from
#'   `formula_par` and `data_par`. There is one model frame
#'   for each nonlinear model parameter listed in
#'   `get_par_names(curve, excess, distr, weekday, link = TRUE)`.
#'   `frame_par[[name]]` retains `terms(formula_par[[name]])` as
#'   an attribute. Model frames correspond rowwise to `endpoints`.
#'   That is, mixed effects data on the fitting window defined
#'   by row `i` of `endpoints` can be found in row `i` of each
#'   model frame.
#' }
#' \item{`frame_append`}{
#'   A data frame preserving the variables from `data_par`
#'   indicated by `append`. Corresponds rowwise to `endpoints`.
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
#' If `debug = TRUE`, then a list containing only
#' `endpoints`, `frame_ts`, `frame_par`, `init`, `tmb_args`, and `call`.
#'
#' @export
#' @importFrom TMB MakeADFun sdreport
#' @useDynLib epigrowthfit
egf <- function(formula_ts,
                formula_par,
                data_ts = parent.frame(),
                data_par = parent.frame(),
                endpoints,
                origin = .Date(0),
                curve = c("logistic", "richards", "exponential", "subexponential", "gompertz"),
                excess = FALSE,
                distr = c("nbinom", "pois"),
                weekday = FALSE,
                method = c("nlminb", "nlm", "Nelder-Mead", "BFGS"),
                na_action = c("fail", "fail"),
                sparse_X = FALSE,
                debug = FALSE,
                init = NULL,
                append = NULL,
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
  na_action <- match.arg(na_action, c("fail", "exclude", "pass"), several.ok = TRUE)
  na_action <- rep_len(na_action, 2L)
  stop_if_not_true_false(sparse_X)
  stop_if_not_true_false(debug)

  frames <- make_frames(
    formula_ts = formula_ts,
    formula_par = formula_par,
    data_ts = data_ts,
    data_par = data_par,
    endpoints = endpoints,
    origin = origin,
    curve = curve,
    excess = excess,
    distr = distr,
    weekday = weekday,
    na_action = na_action,
    init = init,
    append = substitute(append)
  )
  tmb_args <- make_tmb_args(
    frame_ts = frames$frame_ts,
    frame_par = frames$frame_par,
    curve = curve,
    distr = distr,
    excess = excess,
    weekday = weekday,
    sparse_X = sparse_X,
    init = init,
    debug = debug
  )

  if (debug) {
    init_split <- tmb_args$parameters
    Y_init <- attr(init_split, "Y_init")
    if (!has_random(tmb_args$data)) {
      init_split <- init_split["beta"]
    }
    init <- unlist(init_split)
    names(init) <- enum_dupl_string(rep.int(names(init_split), lengths(init_split)))
    out <- list(
      endpoints = frames$endpoints,
      frame_ts = frames$frame_ts,
      frame_par = frames$frame_par,
      init = init,
      tmb_args = tmb_args,
      Y_init = Y_init,
      call = match.call()
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
    endpoints = frames$endpoints,
    frame_ts = frames$frame_ts,
    frame_par = frames$frame_par,
    frame_append = frames$frame_append,
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
