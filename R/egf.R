#' Fit phenomenological models of epidemic growth
#'
#' @description
#' Fits nonlinear mixed effects models of epidemic growth to one
#' or more disease incidence time series.
#'
#' @param model
#'   An \code{"\link{egf_model}"} object specifying a nonlinear model
#'   and a dispersion model to be estimated.
#' @param formula
#'   A \link{formula} of the form \code{x ~ time} or \code{x ~ time | ts}
#'   specifying one or more incidence time series in long format.
#'   \code{time} must evaluate to a \link{numeric} or \link{Date} vector.
#'   Numeric \code{time} is expected to measure time as a number of days
#'   since \emph{the end of} Date \code{origin}. Date \code{time} is
#'   coerced to numeric \code{\link{julian}(time, origin)}. \code{x}
#'   must evaluate to a numeric vector. Within time series, \code{x[i]}
#'   should specify the number of cases observed from \code{time[i-1]}
#'   to \code{time[i]}. Finally, \code{ts} must evaluate to a \link{factor}
#'   such that \code{\link{split}(\link{data.frame}(time, x), ts)}
#'   returns a \link{list} of time series. Formulae \code{x ~ time}
#'   are equivalent to \code{x ~ time | ts} with \code{ts} set equal
#'   to \code{\link{rep}(\link{factor}(1), \link{length}(x))}.
#' @param formula_par
#'   A \link{list} of \link{formula}e of the form \code{par ~ terms},
#'   specifying mixed effects models for nonlinear and dispersion model
#'   parameters using \code{\link[lme4:lmer]{lme4}}-like syntax.
#'   Alternatively, a formula of the form \code{~terms} to be recycled
#'   for all parameters. A list of the parameters for which formulae
#'   can be specified can be retrieved with \code{\link{get_par_names}}.
#'   Specifically, \code{\link{deparse}(par)} must be an element of
#'   \code{\link{get_par_names}(model, link = TRUE)}. The default for
#'   parameters not assigned a formula is \code{~1}.
#' @param data,data_par
#'   \link[=data.frame]{Data frame}s, \link{list}s, or \link{environment}s.
#'   These are searched prior to \link{formula} environments for variables
#'   used in \code{formula} and \code{formula_par}, respectively.
#' @param subset,subset_par
#'   Expressions to be evaluated in \code{data} and \code{data_par},
#'   respectively, provided the latter are data frames. The values,
#'   which must be \link{logical} vectors of suitable length,
#'   index rows of data to be used when fitting (see Details).
#'   The default (\code{\link{NULL}}) is to use all rows.
#' @param na_action
#'   A \link{character} string affecting the handling of \code{\link{NA}}
#'   in \code{x} for \code{formula = x ~ time | ts}.
#'   \code{"fail"} is to throw an error.
#'   \code{"pass"} is to ignore \code{NA} when fitting and replace \code{NA}
#'   when predicting.
#'   Note that \code{NA} in \code{time} and \code{ts} are an error regardless
#'   of \code{na_action}.
#' @param na_action_par
#'   A \link{character} string affecting the handling of \code{\link{NA}}
#'   in \code{formula_par} variables.
#'   \code{"fail"} is to throw an error.
#'   \code{"omit"} is to discard fitting windows with incomplete data.
#' @param endpoints
#'   A \link[=data.frame]{data frame}, \link{list}, or \link{environment}
#'   with variables \code{start} and \code{end}, and any further variables
#'   necessary to evaluate \code{ts} when \code{formula = x ~ time | ts}.
#'   \code{start} and \code{end} must be \link{numeric} or \link{Date}
#'   vectors listing start and end times for all fitting windows. \code{ts}
#'   must evaluate to a \link{factor} indicating the time series in which each
#'   window is found. Within time series, intervals \code{[start[i], end[i]]}
#'   must be disjoint and contain at least two time points from \code{time}.
#' @param control
#'   An \code{"\link{egf_control}"} object specifying control parameters.
#' @param do_fit
#'   A \link{logical} flag. If \code{TRUE}, then \code{egf} returns early
#'   with a \link{list} of optimization inputs.
#' @param se
#'   A \link{logical} flag. If \code{TRUE}, then standard errors on mixed
#'   effects model coefficients are computed and stored for later reuse
#'   by methods.
#' @param init
#'   A \link{numeric} vector to be used as the full parameter vector for
#'   the first likelihood evaluation. The default (\code{\link{NULL}}) is
#'   to accept the internally generated default. Use \code{do_fit = FALSE}
#'   to retrieve this default and other optimization inputs (in particular
#'   \code{Y_init}; see Value), which can often be used to construct a more
#'   informative starting point.
#' @param origin
#'   A \link{Date} specifying a reference time.
#' @param append
#'   An expression indicating variables in \code{data_par} to be preserved
#'   in the returned object for use by methods. It is evaluated similarly
#'   to the \code{select} argument of \code{\link{subset}}. The default
#'   (\code{\link{NULL}}) is to preserve all variables. Currently, usage
#'   is supported only if \code{data_par} is a \link[=data.frame]{data frame}.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' \code{endpoints} and \link[=model.frame]{model frame}s constructed from
#' \code{formula_par} and \code{data_par} are expected to correspond rowwise.
#'
#' Day of week effect estimation (see \code{\link{egf_model}}) requires
#' \code{time} in \code{formula = x ~ time | ts} to evaluate to
#' an \link{integer} or \link{Date} vector with 1-day spacing in
#' all fitting windows.
#'
#' @return
#' If \code{do_fit = TRUE}, then a \link{list} inheriting from \link{class}
#' \code{"egf"}, with elements:
#' \item{endpoints}{
#'   A \link[=data.frame]{data frame} with variables \code{ts}, \code{window},
#'   \code{start}, and \code{end} listing start and end times for all fitting
#'   windows. Rows are ordered by time series and chronologically within time
#'   series. \code{origin} is retained as an \link[=attributes]{attribute}.
#'   Note that \code{start} and \code{end} here are the minimum and maximum of
#'   the set of time points contained in the supplied intervals. Hence supplied
#'   and returned interval endpoints may not match exactly.
#' }
#' \item{frame}{
#'   The time series model frame, constructed from \code{formula} and
#'   \code{data}, with variables \code{ts}, \code{window}, \code{start},
#'   and \code{end}. \code{ts} and \code{window} are \link{factor}s
#'   that can be used to \link{split} \code{frame} by time series and
#'   by fitting window, respectively. Rows are ordered by time series
#'   and chronologically within time series. \code{\link{terms}(formula)}
#'   and \code{origin} are retained as \link{attributes}.
#'   \code{\link{as.integer}(window)} indexes rows of \code{endpoints}.
#'   That is, time series data for the fitting window defined in row
#'   \code{i} of \code{endpoints} can be found in the rows of \code{frame}
#'   for which \code{window = \link{levels}(window)[i]}.
#' }
#' \item{frame_par}{
#'   A \link{list} of mixed effects model frames, constructed from
#'   \code{formula_par} and \code{data_par} using \code{\link{model.frame}}.
#'   There is one model frame for each nonlinear and dispersion model
#'   parameter. \code{frame_par[[name]]} retains
#'   \code{\link{terms}(formula_par[[name]])} as
#'   an \link[=attributes]{attribute}.
#'   Model frames correspond rowwise to \code{endpoints}. That is,
#'   mixed effects data on the fitting window defined in row \code{i}
#'   of \code{endpoints} can be found in row \code{i} of each model frame.
#' }
#' \item{frame_append}{
#'   A \link[=data.frame]{data frame} preserving the variables from
#'   \code{data_par} indicated by \code{append}. Corresponds rowwise
#'   to \code{endpoints}.
#' }
#' \item{model, control}{
#'   Copies of the so-named arguments.
#' }
#' \item{tmb_args}{
#'   A \link{list} of arguments to \code{\link[TMB]{MakeADFun}}
#'   constructed by \code{\link{make_tmb_args}}.
#' }
#' \item{tmb_out}{
#'   The \link{list} output of \code{\link[TMB]{MakeADFun}}
#'   (after optimization).
#' }
#' \item{optimizer_out}{
#'   The \link{list} output of the outer optimizer specified
#'   by \code{control$optimizer}.
#' }
#' \item{init, best}{
#'   Numeric vectors. These are the full parameter vectors
#'   of the first and best likelihood evaluations.
#' }
#' \item{nonrandom}{
#'   An \link{integer} vector indexing components \code{beta} and \code{theta}
#'   of \code{init} and \code{best}.
#' }
#' \item{sdreport}{
#'   If \code{se = TRUE}, then an \code{"sdreport"} object resulting
#'   from \code{\link{sdreport}(tmb_out)}. Otherwise, \code{\link{NULL}}.
#' }
#' \item{nll}{
#'   A \link[=double]{numeric} scalar giving the value of the negative
#'   log Laplace approximation of the marginal likelihood function at
#'   \code{best[nonrandom]}.
#' }
#' \item{call}{
#'   The \code{\link{call}} to \code{egf}, allowing for updates
#'   to the \code{"egf"} object via \code{\link{update}}.
#' }
#' If \code{do_fit = FALSE}, then a \link{list} containing \code{frame},
#' \code{frame_par}, \code{endpoints}, \code{tmb_args}, \code{tmb_out}
#' (\emph{before} optimization), \code{init}, \code{nonrandom}, and
#' \code{call}, as well as a \link{matrix} \code{Y_init} specifying
#' a naive estimate of each nonlinear and dispersion model parameter
#' for each fitting window, corresponding rowwise to \code{endpoints}.
#'
#' @export
#' @useDynLib epigrowthfit
egf <- function(model = egf_model(),
                formula,
                formula_par = list(),
                data = parent.frame(),
                data_par = parent.frame(),
                subset = NULL,
                subset_par = NULL,
                na_action = c("fail", "pass"),
                na_action_par = c("fail", "omit"),
                endpoints,
                control = egf_control(),
                do_fit = TRUE,
                se = FALSE,
                init = NULL,
                origin = .Date(0L),
                append = NULL,
                ...) {
  stop_if_not(
    inherits(model, "egf_model"),
    m = "`model` must inherit from class \"egf_model\". See `?egf_model`."
  )
  stop_if_not(
    inherits(control, "egf_control"),
    m = "`control` must inherit from class \"egf_control\". See `?egf_control`."
  )
  stop_if_not_true_false(do_fit)
  stop_if_not_true_false(se)

  if (control$omp_num_threads > 1L) {
    on <- TMB::openmp(n = NULL)
    if (on > 0L) {
      on.exit(TMB::openmp(n = on))
      TMB::openmp(n = control$omp_num_threads)
    }
  }

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
    origin = origin,
    model = model,
    init = init,
    append = append
  )
  tmb_args <- make_tmb_args(
    frame = frames$frame,
    frame_par = frames$frame_par,
    model = model,
    control = control,
    do_fit = do_fit,
    init = init
  )
  tmb_out <- do.call(TMB::MakeADFun, tmb_args)
  tmb_out$fn <- patch_fn_gr(tmb_out$fn,
    order = 0L,
    inner_optimizer = control$inner_optimizer
  )
  tmb_out$gr <- patch_fn_gr(tmb_out$gr,
    order = 1L,
    inner_optimizer = control$inner_optimizer
  )
  init <- enum_dupl_names(tmb_out$env$par)
  nonrandom <- seq_along(init)
  if (!is.null(tmb_out$env$random)) {
    nonrandom <- nonrandom[-tmb_out$env$random]
  }

  if (!do_fit) {
    out <- frames[c("frame", "frame_par", "endpoints")]
    out <- c(out, list(
      tmb_args = tmb_args,
      tmb_out = tmb_out,
      init = init,
      nonrandom = nonrandom,
      Y_init = attr(tmb_args$parameters, "Y_init"),
      call = match.call()
    ))
    return(out)
  }

  optimizer <- control$optimizer$f
  optimizer_args <- c(
    tmb_out[c("par", "fn", "gr")],
    control$optimizer["control"],
    control$optimizer[["args"]]
  )
  optimizer_out <- do.call(optimizer, optimizer_args)
  best <- enum_dupl_names(tmb_out$env$last.par.best)
  nll <- as.numeric(optimizer_out$value)
  sdreport <- if (se) try(TMB::sdreport(tmb_out))

  out <- c(frames, list(
    model = model,
    control = control,
    tmb_args = tmb_args,
    tmb_out = tmb_out,
    optimizer_out = optimizer_out,
    nonrandom = nonrandom,
    init = init,
    best = best,
    nll = nll,
    sdreport = sdreport,
    call = match.call()
  ))
  class(out) <- c("egf", "list")
  out
}

#' Define an epidemic model
#'
#' Sets flags defining the phenomenological model of epidemic growth
#' to be estimated by \code{\link{egf}}.
#'
#' @param curve
#'   A \link{character} string specifying a model for expected cumulative
#'   incidence as a function of time.
#' @param excess
#'   A \link{logical} flag. If \code{TRUE}, then a constant baseline
#'   mortality rate is estimated. Set to \code{TRUE} if what is observed
#'   is multiple causes mortality rather than disease mortality or disease
#'   incidence.
#' @param family
#'   A \link{character} string specifying a model for observed interval
#'   incidence (i.e., an error distribution).
#' @param day_of_week
#'   An integer flag. If positive, then day of week effects are estimated
#'   as offsets relative to the indicated day of week
#'   (Sunday if \code{day_of_week = 1}, Monday if \code{day_of_week = 2},
#'   and so on). \link[=logical]{Logical} values are coerced to integer.
#'
#' @return
#' A \link{list} inheriting from \link{class} \code{"egf_model"}
#' containing the arguments (after possible matching and coercion).
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
#' Set parameters controlling the behaviour of \code{\link{egf}}.
#'
#' @param optimizer
#'   An \code{"\link{egf_optimizer}"} object, specifying an "outer"
#'   optimization method.
#' @param inner_optimizer
#'   An \code{"\link{egf_inner_optimizer}"} object, specifying an "inner"
#'   optimization method, or a \link{list} of such objects, in which case
#'   the listed methods are tried in order until one succeeds. (If none
#'   succeed, then a warning is issued.)
#' @param trace
#'   An integer flag determining the amount of tracing performed
#'   (see Details). \link[=logical]{Logical} values are coerced to integer.
#' @param profile
#'   A \link{logical} flag. If \code{TRUE}, then fixed effect parameters
#'   are profiled out of the likelihood, which may stabilize optimization
#'   for models with many fixed effects. This feature is experimental,
#'   and in fact may \emph{de}stabilize optimization, as it relies on
#'   assumptions about the optimization problem that are not necessarily
#'   satisfied by the models fit by \code{\link{egf}}.
#' @param sparse_X
#'   A \link{logical} flag. If \code{TRUE}, then the fixed effects
#'   \link[=model.matrix]{design} matrix is constructed in
#'   \link[Matrix:sparseMatrix]{sparse} format.
#' @param omp_num_threads
#'   An integer specifying a number of OpenMP threads to be used
#'   (if supported) when evaluating the objective function.
#'
#' @details
#' \code{trace} affects the amount of information printed during
#' likelihood evaluations:
#' \describe{
#' \item{0}{
#'   Likelihood evaluations are always silent.
#' }
#' \item{1}{
#'   A message is printed whenever a negative log likelihood term
#'   is non-finite or exceeds \code{1e+09}.
#' }
#' \item{3}{
#'   All negative log likelihood terms are printed.
#' }
#' \item{2, 4}{
#'   Equivalent to 1 and 3, but with further printing of the
#'   response matrix \code{Y}. \code{Y[i, j]} is the current
#'   value of nonlinear or dispersion model parameter \code{j}
#'   (link scale) in fitting window \code{i}.
#' }
#' }
#'
#' \code{\link{egf}} passes \code{silent = (trace == 0L)}
#' to \code{\link[TMB]{MakeADFun}}. As a result, nonzero values
#' of \code{trace} have a number of additional side effects:
#' \itemize{
#' \item error messages are printed during function and gradient evaluations;
#' \item the maximum absolute gradient element is printed with each gradient
#' evaluation; and
#' \item trace flags set by \code{\link[TMB]{config}} are turned on.
#' }
#'
#' @return
#' A \link{list} inheriting from \link{class} \code{"egf_control"}
#' containing the arguments (after possible coercion).
#'
#' @export
egf_control <- function(optimizer = egf_optimizer(),
                        inner_optimizer = egf_inner_optimizer(),
                        trace = FALSE,
                        profile = FALSE,
                        sparse_X = FALSE,
                        omp_num_threads = getOption("egf.cores", 1L)) {
  stop_if_not(
    inherits(optimizer, "egf_optimizer"),
    m = "`optimizer` must inherit from class \"egf_optimizer\". See `?egf_optimizer`."
  )
  if (inherits(inner_optimizer, "egf_inner_optimizer")) {
    inner_optimizer <- list(inner_optimizer)
  }
  stop_if_not(
    is.list(inner_optimizer),
    vapply(inner_optimizer, inherits, FALSE, "egf_inner_optimizer"),
    m = paste0(
      "`inner_optimizer` must inherit from class \"egf_inner_optimizer\"\n",
      "or be a list of such objects. See `?egf_inner_optimizer`."
    )
  )
  stop_if_not_true_false(trace, allow_numeric = TRUE)
  trace <- min(4L, max(0L, as.integer(trace))) # coercion to `0:4`
  stop_if_not_true_false(profile)
  stop_if_not_true_false(sparse_X)
  stop_if_not_integer(omp_num_threads, kind = "positive")

  out <- list(
    optimizer = optimizer,
    inner_optimizer = inner_optimizer,
    trace = trace,
    profile = profile,
    sparse_X = sparse_X,
    omp_num_threads = omp_num_threads
  )
  class(out) <- c("egf_control", "list")
  out
}


#' Define an optimization method
#'
#' These two functions link an optimizer with function arguments and
#' control parameters to define an optimization method for use by
#' \code{\link{egf}}. "Outer" and "inner" optimization methods are
#' defined separately by \code{egf_optimizer} and \code{egf_inner_optimizer}.
#'
#' @param f
#'   A \link{function} performing optimization. The outer optimization
#'   permits \code{\link{optim}}, \code{\link{nlminb}}, and \code{\link{nlm}}
#'   and any \code{\link{optim}}-like function. An \code{\link{optim}}-like
#'   function is a function \code{f} such that:
#'   (i) the first three arguments of \code{f} specify an initial parameter
#'   vector, an objective function, and a gradient function, respectively;
#'   (ii) \code{f} accepts \code{control} as a fourth (or later) argument;
#'   and
#'   (iii) \code{f} returns a \link{list} with elements \code{par},
#'   \code{value}, \code{convergence}, and \code{message}.
#'   The inner optimization permits \code{\link{optim}}
#'   and \code{\link[TMB]{newton}} only.
#' @param args
#'   A \link{list} of arguments to \code{f} other than \code{control}.
#'   If \code{f = \link{optim}} and \code{args} does not have \code{method}
#'   as an element, then \code{method = "BFGS"} is appended.
#' @param control
#'   A \link{list} of control parameters to be assigned to argument
#'   \code{control} of \code{f}.
#'
#' @return
#' \code{egf_optimizer} returns a \link{list} inheriting from \link{class}
#' \code{"egf_optimizer"}, with elements:
#' \item{f}{
#'   An \code{\link{optim}}-like \link{function}, typically the result
#'   of wrapping the supplied optimizer to make it \code{\link{optim}}-like.
#' }
#' \item{args}{
#'   The supplied \link{list} of arguments
#'   (after possible deletion of elements with reserved names).
#' }
#' \item{control}{
#'   The supplied \link{list} of control parameters.
#' }
#'
#' \code{egf_inner_optimizer} returns a \link{list} inheriting
#' from \link{class} \code{"egf_inner_optimizer"}, with elements:
#' \item{method}{
#'   A \link{character} string. This is
#'   \code{args$method} if \code{f = \link{optim}} and
#'   \code{"newton"} if \code{f = \link[TMB]{newton}}.
#' }
#' \item{control}{
#'   A \link{list}. This is \code{control} if \code{f = \link{optim}} and
#'   \code{args} (after possible deletion of elements with reserved names)
#'   if \code{f = \link[TMB]{newton}}.
#' }
#'
#' @seealso
#' \code{\link[TMB]{MakeADFun}} for some details about outer and inner
#' optimizations
#'
#' @name egf_optimizer
NULL

#' @rdname egf_optimizer
#' @export
#' @importFrom stats optim nlminb nlm
#' @importFrom TMB newton
egf_optimizer <- function(f = nlminb, args = list(), control = list()) {
  s <- substitute(f)
  optimizer <- f
  stop_if_not(
    is.list(args),
    m = "`args` must be a list."
  )
  stop_if_not(
    is.list(control),
    m = "`control` must be a list."
  )
  if (identical(f, optim)) {
    if (is.null(args$method)) {
      args$method <- "BFGS"
      message(sprintf("`optim` method not specified, using \"%s\".", args$method))
    } else {
      args$method <- match.arg(args$method, eval(formals(optim)$method))
    }
  } else if (identical(f, nlminb)) {
    f <- function(par, fn, gr, control, ...) {
      l <- nlminb(start = par, objective = fn, gradient = gr, control = control, ...)
      l["value"] <- l["objective"]
      l
    }
  } else if (identical(f, nlm)) {
    f <- function(par, fn, gr, control, ...) {
      l <- nlm(f = structure(fn, gradient = gr), p = par, ...)
      l[c("par", "value", "convergence")] <- l[c("estimate", "minimum", "code")]
      l["message"] <- list(NULL)
      l
    }
  } else {
    stop_if_not(
      is.function(f),
      m = sprintf("`%s` must be a function.", s)
    )
    nf <- names(formals(f))
    stop_if_not(
      length(nf) >= 4L,
      nf[1:3] != "...",
      "control" %in% nf[-(1:3)],
      m = sprintf("`formals(%s)` must have configuration outlined in `?egf_optimizer`.", s)
    )
    e <- quote(f(c(1, 1), function(x) sum(x^2), function(x) 2 * x))
    l <- try(eval(e))
    if (inherits(l, "try-error")) {
      stop(
        sprintf("Unable to validate optimizer `%s` because test\n", s),
        sprintf("`%s` produced an error.", sub("^f", s, deparse(e)))
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
    f <- function(par, fn, gr, control, ...) {
      optimizer(par, fn, gr, control = control, ...)
    }
  }
  if (!is.null(names(args))) {
    reserved <- c("par", "fn", "gr", "control", "...", names(formals(optimizer))[1:3])
    args <- args[setdiff(names(args), reserved)]
  }

  out <- list(f = f, args = args, control = control)
  class(out) <- c("egf_optimizer", "list")
  out
}

#' @rdname egf_optimizer
#' @export
#' @importFrom stats optim
#' @importFrom TMB newton
egf_inner_optimizer <- function(f = newton, args = list(), control = list()) {
  stop_if_not(
    is.list(args),
    m = "`args` must be a list."
  )
  stop_if_not(
    is.list(control),
    m = "`control` must be a list."
  )

  if (identical(f, optim)) {
    if (is.null(args$method)) {
      method <- "BFGS"
      message(sprintf("`optim` method not specified, using \"%s\".", method))
    } else {
      method <- match.arg(args$method, eval(formals(optim)$method))
    }
  } else if (identical(f, newton)) {
    method <- "newton"
    if (!is.null(names(args))) {
      reserved <- c("par", "fn", "gr", "he", "env", "...")
      args <- args[setdiff(names(args), reserved)]
    }
    control <- args
  } else {
    stop("`f` is currently restricted to `TMB::newton` and `stats::optim`.")
  }

  out <- list(method = method, control = control)
  class(out) <- c("egf_inner_optimizer", "list")
  out
}
