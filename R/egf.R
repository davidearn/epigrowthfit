#' Fit models of epidemic growth
#'
#' @description
#' Fits nonlinear mixed effects models of epidemic growth to one
#' or more disease incidence time series.
#'
#' @param formula
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
#'   element of `get_par_names(model, link = TRUE)`.
#'   `~1` is the default for parameters not assigned a formula.
#'   Alternatively, `formula_par` may be a formula of the form
#'   `~terms`. In this case, the formula is recycled for all
#'   nonlinear model parameters. Note that "individuals" in each
#'   mixed effects model are fitting windows, and model frames
#'   constructed from `formula_par` and `data_par` are expected
#'   to correspond rowwise to `endpoints` (see Details).
#' @param data,data_par
#'   Data frames, lists, or environments. These are searched prior
#'   to formula environments for variables used in `formula` and
#'   `formula_par`, respectively.
#' @param subset,subset_par
#'   Expressions to be evaluated in `data` and `data_par`, respectively,
#'   specifying observations to be used to fit the mixed effects model.
#'   These must evaluate to logical vectors of suitable length.
#'   Note that `subset` and `subset_par` perform subsetting at
#'   different levels of granularity: each row of `data` corresponds
#'   to a time point in a time series, whereas each row of `data_par`
#'   corresponds to a fitting window (and thus multiple time points)
#'   in a time series.
#' @param na_action
#'   A character string affecting the handling of `NA` in `x`
#'   if `formula = x ~ time | ts`.
#'   `"fail"` is to throw an error.
#'   `"pass"` is to ignore `NA` when fitting and replace `NA`
#'   when predicting.
#'   (TODO:
#'   `"exclude"` is to ignore `NA` when fitting and preserve `NA`
#'   when predicting.)
#'   Note that `NA` in `time` and `ts` are an error regardless
#'   of `na_action`.
#' @param na_action_par
#'   A character string affecting the handling of `NA`
#'   in `formula_par` variables.
#'   `"fail"` is to throw an error.
#'   `"pass"` is to discard fitting windows with incomplete data.
#' @param endpoints
#'   A data frame, list, or environment with variables `start`
#'   and `end`, and any further variables necessary to evaluate
#'   `ts` if `formula = x ~ time | ts`. `start` and `end` must
#'   be numeric or Date vectors listing start and end times
#'   for all fitting windows. `ts` must evaluate to a factor
#'   indicating the time series in which each window is found.
#'   Within time series, intervals `[start[i], end[i]]` must be
#'   disjoint and contain at least two time points from `time`.
#' @param model
#'   An `"egf_model"` object constructed using [egf_model()],
#'   specifying a nonlinear model to be estimated.
#' @param control
#'   An `"egf_control"` object constructed using [egf_control()],
#'   specifying control parameters.
#' @param do_fit
#'   A logical scalar. If `TRUE`, then `egf()` returns early
#'   with a list of optimization inputs.
#' @param se
#'   A logical scalar. If `TRUE`, then standard errors on mixed
#'   effects model coefficients are computed and stored for later
#'   reuse by methods.
#' @param init
#'   A full parameter vector for the first likelihood evaluation.
#'   The default (`NULL`) is to accept the internally generated
#'   default. Use `do_fit = FALSE` to retrieve this default and
#'   other optimization inputs (in particular `Y_init`; see Value),
#'   which can often be used to construct a more informative
#'   starting point.
#' @param origin
#'   A Date specifying a reference time.
#' @param append
#'   An expression indicating variables in `data_par` to be preserved
#'   in the returned `"egf"` object for use by methods.
#'   It is evaluated similarly to the `select` argument of [subset()].
#'   The default (`NULL`) is to preserve all variables.
#'   Currently, usage is supported only if `data_par` is a data frame.
#' @param ...
#'   Unused optional arguments
#'
#' @details
#' If `formula_ts = x ~ time | ts`, then coercion of Date `time` to
#' numeric `julian(time, origin)` assumes that each element `time[i]`
#' can be read as "end of date `time[i]`", so that observation `x[i]`
#' in a given time series is the number of cases observed from the
#' end of date `time[i-1]` to the end of date `time[i]`.
#'
#' To avoid unexpected mismatch between `endpoints` and mixed effects
#' model frames constructed from `formula_par` and `data_par`,
#' keep all necessary `endpoints` variables in `data_par`, and
#' set `endpoints = data_par`.
#'
#' Day of week effect estimation (see [egf_model()]) requires `time`
#' in `formula = x ~ time | ts` to evaluate to an integer
#' (in the sense of `all.equal(time, round(time))`) or Date vector
#' with 1-day spacing in all fitting windows.
#'
#' @return
#' If `do_fit = TRUE`, then a list inheriting from class `"egf"`,
#' with elements:
#' \item{`endpoints`}{
#'   A data frame with variables `ts`, `window`, `start`, and
#'   `end` listing start and end times for all fitting windows.
#'   Rows are ordered by time series and chronologically
#'   within time series. `origin` is retained as an attribute.
#'   Note that `start` and `end` here are the minimum and maximum
#'   of the set of time points contained in the supplied intervals.
#'   Hence supplied and returned interval endpoints may not
#'   match exactly.
#' }
#' \item{`frame_ts`}{
#'   The time series model frame constructed from `formula`
#'   and `data`, with variables `ts`, `window`, `time`, and `x`.
#'   `ts` and `window` are factors that can be used to split
#'   `frame` by time series and by fitting window, respectively.
#'   Rows are ordered by time series and chronologically within
#'   time series. `terms(formula)` and `origin` are retained
#'   as attributes. `as.integer(window)` indexes rows of `endpoints`.
#'   That is, time series data for the fitting window defined in
#'   row `i` of `endpoints` can be found in the rows of `frame_ts`
#'   for which `window = levels(window)[i]`.
#' }
#' \item{`frame_par`}{
#'   A list of mixed effects model frames, constructed from
#'   `formula_par` and `data_par`. There is one model frame
#'   for each nonlinear model parameter listed in
#'   `get_par_names(model, link = TRUE)`. `frame_par[[name]]`
#'   retains `terms(formula_par[[name]])` as an attribute.
#'   Model frames correspond rowwise to `endpoints`. That is,
#'   mixed effects data on the fitting window defined in row `i`
#'   of `endpoints` can be found in row `i` of each model frame.
#' }
#' \item{`frame_append`}{
#'   A data frame preserving the variables from `data_par`
#'   indicated by `append`. Corresponds rowwise to `endpoints`.
#' }
#' \item{`model`, `control`}{
#'   Copies of the so-named arguments.
#' }
#' \item{`tmb_args`}{
#'   A list of arguments to [TMB::MakeADFun()]. See [make_tmb_args()].
#' }
#' \item{`tmb_out`}{
#'   The list output of [TMB::MakeADFun()] after optimization.
#' }
#' \item{`optimizer_out`}{
#'   The list output of the outer optimizer specified by `control$optimizer`.
#' }
#' \item{`init`}{
#'   The full parameter vector of the first likelihood evaluation.
#' }
#' \item{`best`}{
#'   The full parameter vector of the best likelihood evaluation.
#' }
#' \item{`nonrandom`}{
#'   An integer vector indexing the non-random segment of `best`.
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
#'
#' If `do_fit = FALSE`, then a list containing `endpoints`,
#' `frame_ts`, `frame_par`, `init`, `tmb_args`, and `call`
#' (as described above), as well as a matrix `Y_init` specifying
#' a naive estimate of each nonlinear model parameter for each
#' fitting window. `Y_init` and `endpoints` correspond rowwise.
#'
#' @export
#' @importFrom TMB MakeADFun sdreport
#' @useDynLib epigrowthfit
egf <- function(formula,
                formula_par,
                data = parent.frame(),
                data_par = parent.frame(),
                subset = NULL,
                subset_par = NULL,
                na_action = c("fail", "pass"),
                na_action_par = c("fail", "pass"),
                endpoints,
                model = egf_model(),
                control = egf_control(),
                do_fit = TRUE,
                se = FALSE,
                init = NULL,
                origin = .Date(0),
                append = NULL,
                ...) {
  stop_if_not_true_false(do_fit)
  stop_if_not_true_false(se)

  frames <- make_frames(
    formula = formula,
    formula_par = formula_par,
    data = data,
    data_par = data_par,
    subset = substitute(subset),
    subset_par = substitute(subset_par),
    na_action = match.arg(na_action),
    na_action_par = match.arg(na_action_par),
    endpoints = endpoints,
    model = model,
    init = init,
    origin = origin,
    append = append
  )
  tmb_args <- make_tmb_args(
    frame_ts = frames$frame_ts,
    frame_par = frames$frame_par,
    model = model,
    control = control,
    do_fit = do_fit,
    control = control,
    init = init
  )

  if (!do_fit) {
    init_split <- tmb_args$parameters
    Y_init <- attr(init_split, "Y_init")
    if (!has_random(tmb_args$data)) {
      init_split <- init_split["beta"]
    }
    init <- unlist(init_split, use.names = FALSE)
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
  tmb_out$fn <- patch_fn_gr(tmb_out$fn,
    order = 0L,
    inner_optimizer = control$inner_optimizer
  )
  tmb_out$gr <- patch_fn_gr(tmb_out$gr,
    order = 1L,
    inner_optimizer = control$inner_optimizer
  )

  optimizer <- control$optimizer$optimizer
  optimizer_args <- c(
    tmb_out[c("par", "fn", "gr")],
    control$optimizer["control"],
    control$optimizer$args
  )
  optimizer_out <- do.call(optimizer, optimizer_args)

  init <- enum_dupl_names(tmb_out$env$par)
  best <- enum_dupl_names(tmb_out$env$last.par.best)
  nonrandom <- grep("^b\\[", names(best), invert = TRUE)

  nll <- as.numeric(optimizer_out$value)
  nll_func <- function(x = best) as.numeric(tmb_out$fn(x[nonrandom]))
  nll_grad <- function(x = best) as.numeric(tmb_out$gr(x[nonrandom]))

  report <- tmb_out$report(best)
  if (se) {
    sdr <- try(sdreport(tmb_out))
    if (!inherits(sdr, "try-error")) {
      report <- c(report, list(cov = sdr$cov.fixed), split_sdreport(sdr))
    }
  }

  out <- list(
    endpoints = frames$endpoints,
    frame_ts = frames$frame_ts,
    frame_par = frames$frame_par,
    frame_append = frames$frame_append,
    model = model,
    control = control,
    tmb_args = tmb_args,
    tmb_out = tmb_out,
    optimizer_out = optimizer_out,
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

#' Define an epidemic model
#'
#' Sets flags defining the phenomenological model of epidemic growth
#' to be estimated by [egf()].
#'
#' @param curve
#'   A character string specifying a model for expected cumulative
#'   incidence as a function of time.
#' @param excess
#'   A logical scalar. If `TRUE`, then a constant baseline mortality
#'   rate is estimated. Set to `TRUE` if what is observed is multiple
#'   causes mortality rather than disease mortality or disease incidence.
#' @param family
#'   A character string specifying a model for observed interval incidence
#'   (i.e., an error distribution).
#' @param day_of_week
#'   An integer scalar. If `day_of_week > 0`, then day of week effects
#'   are estimated as offsets relative to the indicated day of week
#'   (Sunday if `day_of_week = 1`, Monday if `day_of_week = 2`, and so on).
#'   Logical values of `day_of_week` are coerced to integer.
#'
#' @return
#' A list inheriting from class `"egf_model"` containing the arguments
#' (after possible matching and coercion).
#'
#' @export
egf_model <- function(curve = c("logistic", "richards", "exponential", "subexponential", "gompertz"),
                      excess = FALSE,
                      family = c("nbinom", "pois"),
                      day_of_week = FALSE) {
  curve <- match.arg(curve)
  stop_if_not_true_false(excess)
  family <- match.arg(family)
  stop_if_not_true_false(day_of_week, allow_numeric = TRUE)
  day_of_week <- as.integer(day_of_week)
  day_of_week <- (day_of_week > 0L) * (1L + (day_of_week - 1L) %% 7L) # coercion to `0:7`

  out <- list(
    curve = curve,
    excess = excess,
    family = family,
    day_of_week = day_of_week
  )
  class(out) <- c("egf_model", "list")
  out
}

#' Set control parameters
#'
#' Set parameters controlling the behavior of [egf()].
#'
#' @param optimizer
#'   An `"egf_optimizer"` object returned by [egf_optimizer()],
#'   specifying an "outer" optimization method.
#' @param inner_optimizer
#'   An `"egf_inner_optimizer"` object returned by [egf_inner_optimizer()],
#'   specifying an "inner" optimization method, or a list of such objects,
#'   in which case the listed methods are tried in order until one succeeds,
#'   and a warning is issued if none succeed.
#' @param trace
#'   An integer scalar.
#'   If 0, then likelihood evaluations are always silent.
#'   If 1, then negative log likelihood terms are printed
#'   when they are non-finite or exceed `1e+09`.
#'   If 3, then all negative log likelihood terms are printed.
#'   Values 2 and 4 are equivalent to 1 and 3, respectively,
#'   but with further printing of the response matrix `Y`.
#'   `Y[i, j]` is the current value of nonlinear model
#'   parameter `j` (link scale) in fitting window `i`.
#'   Logical values of `trace` are coerced to integer.
#' @param sparse_X
#'   A logical scalar. If `TRUE`, then the fixed effects
#'   design matrix is constructed in sparse format.
#'
#' @return
#' A list inheriting from class `"egf_control"` containing the arguments
#' (after possible coercion).
#'
#' @export
egf_control <- function(optimizer = egf_optimizer(),
                        inner_optimizer = egf_inner_optimizer(),
                        trace = FALSE,
                        sparse_X = FALSE) {
  stop_if_not(
    inherits(optimizer, "egf_optimizer"),
    m = "`optimizer` must inherit from class \"egf_optimizer\". See `?egf_optimizer`."
  )
  if (inherits(inner_optimizer, "egf_inner_optimizer")) {
    inner_optimizer <- list(inner_optimizer)
  }
  stop_if_not(
    vapply(inner_optimizer, inherits, FALSE, "egf_inner_optimizer"),
    m = paste0(
      "`inner_optimizer` must inherit from class \"egf_inner_optimizer\"\n",
      "or be a list of such objects. See `?egf_inner_optimizer`."
    )
  )
  stop_if_not_true_false(trace, allow_numeric = TRUE)
  trace <- min(4L, max(0L, as.integer(trace))) # coercion to `0:4`
  stop_if_not_true_false(sparse_X)

  out <- list(
    optimizer = optimizer,
    inner_optimizer = inner_optimizer,
    trace = trace,
    sparse_X = sparse_X
  )
  class(out) <- c("egf_control", "list")
  out
}


#' Define an optimization method
#'
#' These two functions link an optimizer with function arguments and
#' control parameters to define an optimization method for use by [egf()].
#' "Outer" and "inner" optimization methods are defined separately
#' by `egf_optimizer()` and `egf_inner_optimizer()` respectively.
#' See [TMB::MakeADFun()] for some details on outer and inner optimizations.
#'
#' @param optimizer
#'   A function performing optimization. The outer optimization permits
#'   [stats::optim()], [stats::nlminb()], [stats::nlm()], and
#'   and any `optim`-like function. An `optim`-like function
#'   is a function `f` such that:
#'   (i) the first three arguments of `f` specify an initial parameter
#'   vector, an objective function, and a gradient function, respectively;
#'   (ii) `f` accepts `control` as a fourth (or later) argument; and
#'   (iii) `f` returns a list with elements `par`, `value`, `convergence`,
#'   and `message`. The inner optimization permits [stats::optim()] and
#'   [TMB::newton()] only.
#' @param args
#'   A list of arguments to `optimizer` other than `control`.
#'   If `optimizer = optim` and `args` does not have `method`
#'   as an element, then `method = "BFGS"` is appended.
#' @param control
#'   A list of control parameters to be assigned to `optimizer` argument
#'   `control`.
#'
#' @details
#'
#'
#'
#' @return
#' `egf_optimizer()` returns a list inheriting
#' from class `"egf_optimizer"`, with elements:
#' \item{`optimizer`}{
#'   An `optim`-like function.
#'   This may by the result of wrapping the supplied optimizer
#'   to make it `optim`-like.
#' }
#' \item{`args`}{
#'   The supplied list of arguments,
#'   (after possible deletion of elements with reserved names).
#' }
#' \item{`control`}{
#'   The supplied list of control parameters.
#' }
#'
#' `egf_inner_optimizer()` returns a list inheriting
#' from class `"egf_inner_optimizer"`, with elements:
#' \item{`method`}{
#'   A character string. This is `args$method` if `optimizer = optim`
#'   and `"newton"` if `optimizer = newton`.
#' }
#' \item{`control`}{
#'   A list. This is `control` if `optimizer = optim` and `args`
#'   (after possible deletion of elements with reserved names)
#'   if `optimizer = newton`.
#' }
#' See handling of `inner.method` and `inner.control`
#' in the body of [TMB::MakeADFun()] for details.
#'
#' @name egf_optimizer
NULL

#' @rdname egf_optimizer
#' @export
#' @importFrom stats optim nlminb nlm
#' @importFrom TMB newton
egf_optimizer <- function(optimizer = nlminb, args = list(), control = list()) {
  s <- substitute(optimizer)
  o <- optimizer
  stop_if_not(
    is.list(args),
    m = "`args` must be a list."
  )
  stop_if_not(
    is.list(control),
    m = "`control` must be a list."
  )
  if (identical(optimizer, optim)) {
    if (is.null(args$method)) {
      args$method <- "BFGS"
      message(sprintf("`optim` method not specified, using \"%s\".", method))
    } else {
      args$method <- match.arg(args$method, eval(formals(optim)$method))
    }
  }
  if (identical(optimizer, nlminb)) {
    optimizer <- function(par, fn, gr, control, ...) {
      l <- nlminb(start = par, objective = fn, gradient = gr, control = control, ...)
      l["value"] <- l["objective"]
      l
    }
  } else if (identical(optimizer, nlm)) {
    optimizer <- function(par, fn, gr, control, ...) {
      l <- nlm(f = structure(fn, gradient = gr), p = par, ...)
      l[c("par", "value", "convergence")] <- l[c("estimate", "minimum", "code")]
      l["message"] <- list(NULL)
      l
    }
  } else {
    stop_if_not(
      is.function(optimizer),
      m = sprintf("`%s` must be a function.", s)
    )
    nf <- names(formals(optimizer))
    stop_if_not(
      length(nf) >= 4L,
      nf[1:3] != "...",
      "control" %in% nf[-(1:3)],
      m = sprintf("`formals(%s)` must have configuration outlined in `?egf_optimizer.`", s)
    )
    e <- quote(optimizer(c(1, 1), function(x) sum(x^2), function(x) 2 * x))
    l <- try(eval(e))
    if (inherits(l, "try-error")) {
      stop(
        sprintf("Unable to validate optimizer `%s` because test\n", s),
        sprintf("`%s` produced an error.", sub("^optimizer", s, deparse(e)))
      )
    }
    stop_if_not(
      is.list(l),
      (nl <- c("par", "value", "convergence", "message")) %in% names(l),
      m = paste0(
        sprintf("`%s` must return a list with elements\n", s),
        paste(sprintf("`%s`", nl), collapse = ", "), ",\n",
        sprintf("but _does not_ for test `%s`.", sub("^f", s, deparse(e)))
      )
    )
    optimizer <- function(par, fn, gr, control, ...) {
      o(par, fn, gr, control = control, ...)
    }
  }
  if (!is.null(names(args))) {
    reserved <- c("par", "fn", "gr", "control", "...", names(formals(o))[1:3])
    args <- args[setdiff(names(args), reserved)]
  }

  out <- list(optimizer = optimizer, args = args, control = control)
  class(out) <- c("egf_optimizer", "list")
  out
}

#' @rdname egf_optimizer
#' @export
#' @importFrom stats optim
#' @importFrom TMB newton
egf_inner_optimizer <- function(optimizer = newton, args = list(), control = list()) {
  stop_if_not(
    is.list(args),
    m = "`args` must be a list."
  )
  stop_if_not(
    is.list(control),
    m = "`control` must be a list."
  )

  if (identical(optimizer, optim)) {
    if (is.null(args$method)) {
      method <- "BFGS"
      message(sprintf("`optim` method not specified, using \"%s\".", method))
    } else {
      method <- match.arg(args$method, eval(formals(optim)$method))
    }
  } else if (identical(optimizer, newton)) {
    method <- "newton"
    if (!is.null(names(args))) {
      reserved <- c("par", "fn", "gr", "he", "env", "...")
      args <- args[setdiff(names(args), reserved)]
    }
    control <- args
  } else {
    stop("`optimizer` is currently restricted to `TMB::newton` and `stats::optim`.")
  }

  out <- list(method = method, control = control)
  class(out) <- c("egf_inner_optimizer", "list")
  out
}
