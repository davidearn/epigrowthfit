#' Fit phenomenological models of epidemic growth
#'
#' @description
#' Fits nonlinear mixed effects models of epidemic growth to one
#' or more disease incidence time series.
#'
#' @param object
#'   An \R object for which an \code{egf} method exists,
#'   typically an \code{"\link{egf_model}"} object specifying
#'   a nonlinear model and a dispersion model to be estimated.
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
#'   for all parameters. A list of parameters for which formulae may
#'   be specified can be retrieved with \code{\link{get_par_names}}.
#'   Specifically, \code{\link{deparse}(par)} must be an element of
#'   \code{\link{get_par_names}(object, link = TRUE)}. The default for
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
#' @param priors
#'   An \link{list} of \link{formula}e of the form \code{par ~ f(...)},
#'   specifying priors on nonlinear and dispersion model parameters.
#'   \code{\link{deparse}(par)} must be an element of
#'   \code{\link{get_par_names}(object, link = TRUE)}.
#'   \code{f(...)} must be a \link{call} to a \link[=egf_prior]{prior function}
#'   with arguments specifying suitable hyperparameters.
#'   This call is evaluated in the corresponding formula environment.
#'   (Currently, the only implemented prior is Gaussian. As a result,
#'   formulae must have the form \code{par ~ \link{Normal}(mu, sigma)}.)
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
#' @param map
#'   A \link{factor} of length \code{\link{length}(init)}. Elements
#'   of the full parameter vector corresponding to \code{\link{NA}}
#'   in \code{map} are fixed at their initial values, rather than
#'   estimated. Elements corresponding to a common factor level
#'   are constrained to have a common value during estimation.
#'   \link[=logical]{Logical} values of \code{map} are accepted
#'   and coerced to factor internally with
#'   \code{map <- \link{replace}(\link{gl}(\link{length}(map), 1L), map, NA)}.
#' @param origin
#'   A \link{Date} specifying a reference time.
#' @param append
#'   An expression indicating variables in \code{data_par} to be preserved
#'   in the returned object for use by methods. It is evaluated similarly
#'   to the \code{select} argument of \code{\link{subset}}. The default
#'   (\code{\link{NULL}}) is to preserve all variables. Currently, usage
#'   is supported only if \code{data_par} is a \link[=data.frame]{data frame}.
#' @param ...
#'   Arguments passed to methods by the generic function.
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
#' \item{model}{
#'   A copy of argument \code{object}.
#' }
#' \item{priors}{
#'   A \link{list} of \code{"\link{egf_prior}"} objects,
#'   obtained by evaluating the right hand side of each
#'   \link{formula} listed in the so-named argument.
#' }
#' \item{control}{
#'   A copy of the so-named argument.
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
#'   An \link{integer} vector indexing the elements of \code{init}
#'   and \code{best} that are not random effects.
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
#'   The \link{call} to \code{egf}, allowing for updates to the
#'   \code{"egf"} object via \code{\link{update}}.
#' }
#' If \code{do_fit = FALSE}, then a \link{list} inheriting from
#' \link{class} \code{"egf_no_fit"} containing \code{model},
#' \code{endpoints}, \code{frame}, \code{frame_par}, \code{tmb_args},
#' \code{tmb_out} (\emph{before} optimization), \code{init},
#' \code{nonrandom}, and \code{call}, as well as a \link[=double]{numeric}
#' \link{matrix} \code{Y_init} specifying a naive estimate
#' of each nonlinear and dispersion model parameter for each
#' fitting window, corresponding rowwise to \code{endpoints}.
#'
#' @export
#' @useDynLib epigrowthfit
egf <- function(object, ...) {
  UseMethod("egf", object)
}

#' @rdname egf
#' @export
egf.egf_model <- function(object,
                          formula,
                          formula_par = list(),
                          data = parent.frame(),
                          data_par = parent.frame(),
                          subset = NULL,
                          subset_par = NULL,
                          na_action = c("fail", "pass"),
                          na_action_par = c("fail", "omit"),
                          endpoints,
                          priors = list(),
                          control = egf_control(),
                          do_fit = TRUE,
                          se = FALSE,
                          init = NULL,
                          map = NULL,
                          origin = .Date(0L),
                          append = NULL,
                          ...) {
  stop_if_not(
    inherits(control, "egf_control"),
    m = "`control` must inherit from class \"egf_control\". See `?egf_control`."
  )
  stop_if_not_true_false(do_fit)
  stop_if_not_true_false(se)
  stop_if_not_Date(origin)

  if (control$omp_num_threads > 0L) {
    on <- TMB::openmp(n = NULL)
    if (on > 0L) {
      on.exit(TMB::openmp(n = on))
      TMB::openmp(n = control$omp_num_threads)
    }
  }

  frames <- make_frames(
    model = object,
    formula = formula,
    formula_par = formula_par,
    data = data,
    data_par = data_par,
    subset = substitute(subset),
    subset_par = substitute(subset_par),
    na_action = match.arg(na_action),
    na_action_par = match.arg(na_action_par),
    endpoints = endpoints,
    init = init,
    origin = origin,
    append = append
  )
  priors <- make_priors(priors = priors, model = object)
  tmb_args <- make_tmb_args(
    model = object,
    frame = frames$frame,
    frame_par = frames$frame_par,
    priors = priors,
    control = control,
    do_fit = do_fit,
    init = init,
    map = map
  )
  tmb_out <- do.call(TMB::MakeADFun, tmb_args)
  tmb_out$fn <- patch_fn(tmb_out$fn, inner_optimizer = control$inner_optimizer)
  tmb_out$gr <- patch_gr(tmb_out$gr, inner_optimizer = control$inner_optimizer)
  init <- enum_dupl_names(tmb_out$env$par)
  nonrandom <- seq_along(init)
  if (!is.null(tmb_out$env$random)) {
    nonrandom <- nonrandom[-tmb_out$env$random]
  }

  if (!do_fit) {
    out <- frames[c("endpoints", "frame", "frame_par")]
    out <- c(out, list(
      model = object,
      tmb_args = tmb_args,
      tmb_out = tmb_out,
      init = init,
      nonrandom = nonrandom,
      Y_init = attr(tmb_args$parameters, "Y_init"),
      call = match.call()
    ))
    class(out) <- c("egf_no_fit", "list")
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
    model = object,
    priors = priors,
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
  stop_if_not_integer(omp_num_threads, "positive")

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

#' Define a parallelization method
#'
#' Defines instructions for parallelization by linking a method with options.
#'
#' @param method
#'   A \link{character} string indicating a method of parallelization.
#'   \code{"\link[=lapply]{serial}"} indicates no parallelization.
#'   \code{"\link[parallel:mclapply]{multicore}"} indicates forking.
#'   It is intended for use from a terminal rather than a GUI.
#'   On Windows, it is equivalent to \code{"serial"}.
#'   \code{"\link[parallel:parLapply]{snow}"} indicates socket clusters.
#'   It supports parallelization on both Unix-alikes and Windows.
#' @param outfile
#'   A \link{character} string indicating a file path where
#'   console output should be diverted. \code{\link{NULL}} is
#'   equivalent to \code{\link{nullfile}()}. \code{""} means
#'   no diversion. If \code{method = "snow"}, then diversion
#'   may be necessary to view output.
#' @param cores
#'   A positive integer indicating a number of worker processes
#'   to spawn when \code{parallel != "serial"}. The maximum is
#'   typically \code{\link[parallel]{detectCores}(TRUE, FALSE)}.
#' @param options
#'   A \link{list} of optional arguments to
#'   \code{\link[parallel]{mclapply}} (\code{method = "multicore"}) or
#'   \code{\link[parallel]{makePSOCKcluster}} (\code{method = "snow"}).
#' @param cl
#'   An optional \link[parallel:makePSOCKcluster]{socket cluster}
#'   to be used when \code{parallel = "snow"}. If non-\code{\link{NULL}},
#'   then \code{outfile}, \code{cores}, and \code{options} are ignored.
#'
#' @return
#' A \link{list} inheriting from \link{class} \code{"egf_parallel"}
#' containing the arguments (after possible matching and substitution).
#'
#' @export
egf_parallel <- function(method = c("serial", "multicore", "snow"),
                         outfile = "",
                         cores = getOption("egf.cores", 1L),
                         options = list(),
                         cl = NULL) {
  method <- match.arg(method)
  if (is.null(outfile)) {
    outfile <- nullfile()
  } else {
    stop_if_not_string(outfile)
  }
  if (method == "serial") {
    cores <- 1L
    options <- list()
    cl <- NULL
  } else if (method == "multicore" || (method == "snow" && is.null(cl))) {
    stop_if_not_integer(cores, "positive")
    stop_if_not(is.list(options),
      m = "`options` must be a list."
    )
    if (method == "multicore") {
      options$mc.cores <- cores
    } else {
      options$names <- cores
      options$outfile <- outfile
    }
  } else {
    stop_if_not(
      inherits(cl, "SOCKcluster"),
      m = "`cl` must inherit from class \"SOCKcluster\". See `?parallel::makePSOCKcluster`."
    )
    cores <- length(cl)
    options <- list()
  }

  out <- list(
    method = method,
    outfile = outfile,
    cores = cores,
    options = options,
    cl = cl
  )
  class(out) <- c("egf_parallel", "list")
  out
}

#' Simulate incidence time series for tests
#'
#' \code{egf_simulate} simulates daily incidence time series following
#' specified nonlinear and dispersion models, with parameters (optionally)
#' varying between time series according to a random intercept model
#' \code{~(1 | ts)}. \code{egf.egf_simulate} estimates the generative
#' model from the simulated time series.
#'
#' @param N
#'   A positive integer indicating a number of time series.
#' @param model
#'   An \code{"\link{egf_model}"} object specifying a nonlinear model
#'   and a dispersion model to be simulated.
#' @param mu
#'   A \link{numeric} vector listing mean nonlinear and dispersion model
#'   parameter values (link scale). It is assumed that elements are ordered
#'   as in \code{\link{get_par_names}(model, link = TRUE)}.
#' @param Sigma
#'   A symmetric positive semidefinite \link{numeric} \link{matrix},
#'   to be used as the covariance matrix corresponding to \code{mu}.
#'   The default (\code{\link{NULL}}) is equivalent to a zero matrix.
#' @param tol
#'   A non-negative number indicating a tolerance for lack of positive
#'   semidefiniteness of \code{Sigma} (see \code{\link[MASS]{mvrnorm}}).
#' @param tmax
#'   A positive number. Time series generated from a model of
#'   cumulative incidence without an inflection point span this many days.
#' @param cstart
#'   A number indicating a threshold value of cumulative incidence.
#'   Left endpoints of suggested fitting windows are those times when
#'   cumulative incidence first exceeds this threshold. (If it is
#'   never exceeded, then the first time point is used, with a warning.)
#' @param origin
#'   A \link{Date} specifying a reference time.
#' @param object
#'   An \code{"egf_simulate"} object supplying simulated time series
#'   and specifying a generative model to be estimated.
#' @param ...
#'   Optional arguments passed to \code{\link{egf.egf_model}},
#'   such as \code{priors} and \code{control}.
#'
#' @details
#' For models of cumulative incidence with an inflection point
#' (\code{model$curve} equal to \code{"logistic"}, \code{"richards"},
#' or \code{"gompertz"}), simulated time series run from 0 days
#' to \code{ceiling(tinfl)+1} days, where \code{tinfl} is the
#' inflection time.
#'
#' For models of cumulative incidence \emph{without} an inflection point
#' (\code{model$curve} equal to \code{"exponential"}
#' or \code{"subexponential"}), simulated time series run from 0 days
#' to \code{tmax} days.
#'
#' Suggested fitting windows start at a time determined by \code{cstart}
#' and end at the last time point.
#'
#' @return
#' \code{egf_simulate} returns a \link{list} inheriting from \link{class}
#' \code{"egf_simulate"}, with elements:
#' \item{model, mu, Sigma, origin}{
#'   Copies of the so-named arguments.
#' }
#' \item{Y}{
#'   A \link[=double]{numeric} \link{matrix} with \code{N} rows and
#'   \code{length(mu)} columns listing the nonlinear and dispersion
#'   model parameter values underlying each time series.
#'   If \code{Sigma} is \code{\link{NULL}}, then the rows of \code{Y}
#'   are all \code{mu}.
#'   If \code{Sigma} is non-\code{\link{NULL}}, then \code{Y} is the
#'   result of \code{\link[MASS]{mvrnorm}(N, mu, Sigma, tol)}.
#' }
#' \item{formula}{
#'   A \link{formula} expressing how to locate the simulated time series
#'   in \code{data}. The formula is always \code{x ~ time | ts}.
#' }
#' \item{formula_par}{
#'   A \link{formula} specifying the generative model.
#'   If \code{Sigma = \link{NULL}}, then the formula is \code{~1}
#'   if \code{N = 1} and \code{~ts} if \code{N > 1}.
#'   Otherwise, the formula is \code{~(1 | ts)}.
#' }
#' \item{data}{
#'   A \link[=data.frame]{data frame} with variables \code{ts}, \code{time},
#'   and \code{x} storing \code{N} simulated time series in long format.
#'   \code{ts} is a \link{factor} with \code{N} \link{levels} grouping
#'   the rows by time series. \code{time} is a \link[=double]{numeric}
#'   vector listing time points as numbers of days since \emph{the end of}
#'   Date \code{origin}. \code{x} is an \link{integer} vector such that,
#'   within time series, \code{x[i]} is the number of cases observed between
#'   \code{time[i-1]} and \code{time[i]}.
#' }
#' \item{data_par}{
#'   A \link[=data.frame]{data frame} with \code{N} rows defining variables
#'   used in \code{formula_par}. (Hence it is either empty or has one variable
#'   \code{ts}.)
#' }
#' \item{endpoints}{
#'   A \link[=data.frame]{data frame} with \code{N} rows and variables
#'   \code{ts}, \code{start}, and \code{end} suggesting a fitting window
#'   for each simulated time series.
#' }
#' \item{actual}{
#'   A \link[=double]{numeric} vector giving the full parameter vector
#'   corresponding to the generative model. When estimating this model
#'   from \code{data}, \code{\link{egf}} output should be compared against
#'   \code{actual}. More precisely, if \code{m} is the \code{"\link{egf}"}
#'   object, then \code{m$best} estimates \code{actual}.
#' }
#' \item{nonrandom}{
#'   An \link{integer} vector indexing the elements of \code{actual}
#'   that are not random effects.
#' }
#' \item{call}{
#'   The \link{call} to \code{egf_simulate}, allowing for updates
#'   to the \code{"egf_simulate"} object via \code{\link{update}}.
#' }
#'
#' @examples
#' model <- egf_model(curve = "logistic", family = "nbinom")
#'
#' r <- log(2) / 20
#' tinfl <- 160
#' K <- 25000
#' nbdisp <- 50
#'
#' mu <- log(c(r, tinfl, K, nbdisp))
#' Sigma <- diag(rep_len(0.25, length(mu)))
#'
#' set.seed(202737L)
#' sim <- egf_simulate(N = 20L, model = model, mu = mu, Sigma = Sigma)
#' fit <- egf(sim)
#'
#' pp <- data.frame(actual = sim$actual, fitted = fit$best)
#' pp[fit$nonrandom, ]
#'
#' @name egf_simulate
NULL

#' @rdname egf_simulate
#' @export
#' @importFrom stats rpois rnbinom
egf_simulate <- function(N = 1L, model, mu, Sigma = NULL, tol = 1e-06,
                         tmax = 100, cstart = 0, origin = .Date(0L)) {
  stop_if_not(
    requireNamespace("MASS", quietly = TRUE),
    m = wrap(
      "`MASS::mvrnorm` is needed, but `MASS` is not installed. ",
      "Install it by running `install.packages(\"MASS\")`, then try again."
    )
  )
  stop_if_not_integer(N, "positive")
  stop_if_not(
    inherits(model, "egf_model"),
    m = "`model` must inherit from class \"egf_model\". See `?egf_model`."
  )
  par_names <- get_par_names(model, link = TRUE)
  p <- length(par_names)
  stop_if_not(
    is.numeric(mu),
    length(mu) == p,
    is.finite(mu),
    m = sprintf("`mu` must be a finite numeric vector of length %d.", p)
  )
  names(mu) <- par_names
  if (!is.null(Sigma)) {
    stop_if_not(
      is.matrix(Sigma),
      is.numeric(Sigma),
      dim(Sigma) == p,
      is.finite(Sigma),
      isSymmetric(Sigma),
      m = "`Sigma` must be a symmetric finite numeric matrix with `length(mu)` rows."
    )
    dimnames(Sigma) <- rep_len(list(par_names), 2L)
  }
  stop_if_not_number(tol, "nonnegative")
  stop_if_not_number(tmax, "positive")
  stop_if_not_number(cstart)
  stop_if_not_Date(origin)

  if (is.null(Sigma)) {
    Y <- matrix(mu, nrow = N, ncol = p, byrow = TRUE, dimnames = list(NULL, par_names))
  } else {
    Y <- MASS::mvrnorm(N, mu = mu, Sigma = Sigma, tol = tol)
  }

  ## Function simulating observations given time points and
  ## nonlinear and dispersion model parameter values
  do_simulate <- function(time, par) {
    ## Log predicted cumulative incidence
    log_curve <- switch(model$curve,
      exponential = {
        r <- exp(par[["log(r)"]])
        log_c0 <- par[["log(c0)"]]
        log_c0 + r * time
      },
      subexponential = {
        alpha <- exp(par[["log(alpha)"]])
        c0 <- exp(par[["log(c0)"]])
        one_minus_p <- 1 - plogis(par[["logit(p)"]])
        log(c0^one_minus_p + one_minus_p * alpha * time) / one_minus_p
      },
      gompertz = {
        alpha <- exp(par[["log(alpha)"]])
        log_c0 <- par[["log(c0)"]]
        log_K <- par[["log(K)"]]
        log_K + exp(-alpha * time) * (log_c0 - log_K)
      },
      logistic = {
        r <- exp(par[["log(r)"]])
        tinfl <- exp(par[["log(tinfl)"]])
        log_K <- par[["log(K)"]]
        log_K - log1p(exp(-r * (time - tinfl)))
      },
      richards = {
        r <- exp(par[["log(r)"]])
        tinfl <- exp(par[["log(tinfl)"]])
        log_K <- par[["log(K)"]]
        a <- exp(par[["log(a)"]])
        log_K - log1p(a * exp(-r * a * (time - tinfl))) / a
      }
    )
    ## Predicted interval incidence
    cases <- diff(exp(log_curve))
    if (model$excess) {
      ## Incrementing with baseline incidence
      b <- exp(par[["log(b)"]])
      cases <- cases + b * diff(time)
    }
    if (model$day_of_week > 0L) {
      ## Cyclic scaling to model day of week effects
      Date1 <- origin + 1L + time[1L]
      day1 <- julian(Date1, origin = .Date(2L + model$day_of_week)) %% 7L
      w <- c(1, exp(par[sprintf("log(w%d)", 1:6)]))
      cases <- cases * rep_len(w[1L + (day1 + 0:6) %% 7L], length(cases))
    }
    ## Simulated interval incidence
    x <- switch(model$family,
      pois = rpois(length(cases), lambda = cases),
      nbinom = rnbinom(length(cases), mu = cases, size = exp(par[["log(nbdisp)"]]))
    )
    c(NA, x)
  }

  ## Time points ending one day after the inflection time in the
  ## cumulative incidence curve
  ## (or, in the absence of an inflection point, at time `tmax`)
  if (model$curve %in% c("exponential", "subexponential")) {
    time <- rep_len(list(seq.int(from = 0, to = tmax, by = 1)), N)
  } else {
    if (model$curve == "gompertz") {
      tinfl <- log(Y[, "log(K)"] - Y[, "log(c0)"]) / exp(Y[, "log(alpha)"])
    } else {
      tinfl <- exp(Y[, "log(tinfl)"])
    }
    time <- Map(seq.int, from = 0, to = ceiling(tinfl) + 1)
  }
  ## Simulated observations
  x <- Map(do_simulate, time = time, par = lapply(seq_len(N), function(i) Y[i, ]))

  formula <- x ~ time | ts
  data <- data.frame(
    ts = rep.int(gl(N, 1L), lengths(time)),
    time = unlist(time, FALSE, FALSE),
    x = unlist(x, FALSE, FALSE)
  )

  f <- function(x) {
    l <- c(0L, cumsum(x[-1L])) > cstart
    if (any(l)) which.max(l) else NA_integer_
  }
  endpoints <- data.frame(
    ts = gl(N, 1L),
    start = mapply(`[`, time, vapply(x, f, 0L)),
    end = vapply(time, max, 0)
  )
  if (anyNA(endpoints$start)) {
    argna <- which(is.na(endpoints$start))
    warning(
      wrap("Threshold `cstart` not exceeded in these time series:"), "\n\n",
      paste(sprintf("  %*d", nchar(N), argna), collapse = "\n"), "\n\n",
      wrap("Corresponding fitting windows contain all time points (for better or for worse).")
    )
    endpoints$start[argna] <- 0
  }

  if (is.null(Sigma)) {
    if (N == 1L) {
      formula_par <- ~1
      data_par <- endpoints[integer(0L)]
      actual <- mu
      names(actual) <- enum_dupl_string(rep_len("beta", p))
    } else {
      formula_par <- ~ts
      data_par <- endpoints["ts"]
      actual <- rep_len(0, p * N)
      actual[seq.int(from = 1L, by = N, length.out = p)] <- mu
      names(actual) <- enum_dupl_string(rep_len("beta", p))
    }
  } else {
    formula_par <- ~(1 | ts)
    data_par <- endpoints["ts"]
    R <- chol(cov2cor(Sigma))
    iR <- upper.tri(R, diag = TRUE)
    R[iR] <- R[iR] * rep.int(1 / diag(R), seq_len(p))
    l <- list(
      beta = mu,
      b = t(Y) - mu,
      theta = c(0.5 * log(diag(Sigma)), R[upper.tri(R, diag = FALSE)])
    )
    actual <- unlist(l, FALSE, FALSE)
    names(actual) <- enum_dupl_string(rep.int(names(l), lengths(l)))
  }
  environment(formula) <- environment(formula_par) <- .GlobalEnv

  out <- list(
    model = model,
    mu = mu,
    Sigma = Sigma,
    origin = origin,
    Y = Y,
    formula = formula,
    formula_par = formula_par,
    data = data,
    data_par = data_par,
    endpoints = endpoints,
    actual = actual,
    nonrandom = grep("^(beta|theta)\\[", names(actual)),
    call = match.call()
  )
  class(out) <- c("egf_simulate", "list")
  out
}

#' @rdname egf_simulate
#' @export
egf.egf_simulate <- function(object, ...) {
  object$object <- object$model
  s <- c("object", "formula", "formula_par", "data", "data_par", "endpoints", "origin")
  args <- object[s]
  dots <- list(...)
  if (length(dots) > 0L && !is.null(nd <- names(dots))) {
    args <- c(args, dots[match(nd, s, 0L) == 0L])
  }
  do.call(egf, args)
}
