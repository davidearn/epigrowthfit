#' Fit nonlinear mixed effects models of epidemic growth
#'
#' @description
#' Fits nonlinear mixed effects models of epidemic growth
#' to collections of one or more disease incidence time series.
#'
#' @param model
#'   An \code{"\link{egf_model}"} object specifying a top level nonlinear model
#'   to be estimated.
#' @param formula
#'   A \link{formula} of the form \code{cbind(time, x) ~ ts}
#'   specifying one or more disease incidence time series in long format.
#'   \code{ts} must evaluate to a \link{factor}
#'   (insofar that \code{\link{as.factor}(ts)} is a factor)
#'   grouping the data by time series.
#'   \code{time} must evaluate to a \link{numeric} vector
#'   that is increasing within \link{levels} of \code{ts}.
#'   \link{Date} and \link{POSIXt} vectors are tolerated and
#'   coerced to numeric with \code{\link{julian}(time)}.
#'   Finally, \code{x} must evaluate to a non-negative numeric vector
#'   with \code{x[i]} equal to the number of cases observed over the
#'   interval \code{(time[i-1], time[i]]}.
#'   Edge cases like \code{x[1]} are ignored internally.
#'   Nonintegral elements of \code{x} are rounded with a warning.
#'   \code{formula = cbind(time, x) ~ 1} can be supplied
#'   when there is only one time series; it is equivalent
#'   to \code{formula = cbind(time, x) ~ ts} with \code{ts}
#'   evaluating to \code{\link{rep}(factor(1), length(x))}.
#' @param formula_windows
#'   A \link{formula} of the form \code{cbind(start, end) ~ ts_windows}
#'   specifying disjoint fitting windows \code{(start, end]} in long format.
#'   Where \code{formula = cbind(time, x) ~ ts}, observation \code{x[i]}
#'   is associated with window \code{(start[j], end[j]]} if and only if
#'   \code{time[i-1] >= start[j]},
#'   \code{time[i] <= end[j]}, and
#'   \code{ts[i] == ts_windows[j]}.
#' @param formula_parameters
#'   A \link{list} of \link{formula}e of the form \code{parameter ~ terms}
#'   specifying mixed effects models for top level nonlinear model parameters,
#'   using \code{\link[lme4:lmer]{lme4}}-like syntax.
#'   Alternatively, a formula of the form \code{~terms} to be recycled for
#'   all parameters. A list of parameters for which formulae may be specified
#'   can be retrieved with \code{\link{egf_get_names_top}}.
#'   Specifically, \code{\link{deparse}(parameter)} must be an element of
#'   \code{\link{egf_get_names_top}(model, link = TRUE)}.
#'   The default for parameters not assigned a formula is \code{~1}.
#' @param formula_priors_top,formula_priors_bottom
#'   \link[=list]{List}s of \link{formula}e of the form \code{parameter ~ prior}
#'   defining priors on top level nonlinear model parameters or bottom level
#'   mixed effects model parameters, i.e., elements of segments \code{beta}
#'   and \code{theta} of the full parameter vector \code{c(beta, theta, b)}.
#'   \code{prior} must be a \link{call}
#'   to a valid \link[=egf_prior]{prior function}
#'   with arguments specifying suitable hyperparameters
#'   (e.g., \code{\link{Normal}(mu = 0, sigma = 1)}).
#'   In the top level case, \code{deparse(parameter)} must be an element
#'   of \code{\link{egf_get_names_top}(model, link = TRUE)},
#'   and hyperparameters supplied on the right hand side must have length 1.
#'   In the bottom level case, \code{parameter} must be \code{beta}
#'   or \code{theta} or a call to \code{\link{[}} or \code{\link{[[}}
#'   subsetting \code{beta} or \code{theta} (e.g., \code{beta[index]},
#'   where \code{index} is any valid index vector for \code{beta}),
#'   and hyperparameters are recycled to the length of the indicated subset.
#'   All expressions \code{prior} and \code{index} are evaluated in the
#'   corresponding formula environment.
#' @param data,data_windows
#'   \link[=data.frame]{Data frame}s, \link{list}s, or \link{environment}s
#'   to be searched for variables named in the corresponding \link{formula}e.
#'   (\code{formula_parameters} uses \code{data_windows}.)
#'   Formula environments are searched for variables not found here.
#' @param subset,subset_windows
#'   Expressions to be evaluated in the corresponding data frame.
#'   Formula environments are searched for variables not found there.
#'   The result should be a valid index vector for the rows of the
#'   data frame (see \code{\link{[.data.frame}}).
#'   Rows that are indexed are processed further; rows that are not
#'   are discarded.
#'   The default (\code{\link{NULL}}) is to preserve all rows.
#' @param na_action
#'   A \link{character} string affecting the handling of \code{\link{NA}}
#'   in \code{x} if \code{formula = cbind(time, x) ~ ts}.
#'   \code{"fail"} is to throw an error.
#'   \code{"pass"} is to ignore \code{NA} when fitting and replace \code{NA}
#'   when predicting.
#'   Note that \code{NA} in \code{time} and \code{ts} are always an error.
#' @param na_action_windows
#'   A \link{character} string affecting the handling of \code{\link{NA}}
#'   in \code{formula_windows} and \code{formula_parameters} variables.
#'   \code{"fail"} is to throw an error.
#'   \code{"omit"} is to discard incomplete rows of data.
#' @param control
#'   An \code{"\link{egf_control}"} object specifying control parameters.
#' @param fit
#'   A \link{logical} flag. If \code{TRUE}, then \code{egf} returns early
#'   with a \link{list} of optimization inputs.
#' @param se
#'   A \link{logical} flag. If \code{TRUE}, then standard errors on mixed
#'   effects model coefficients are computed and stored for later reuse
#'   by methods.
#' @param init
#'   A \link{numeric} vector to be used as the full parameter vector
#'   \code{c(beta, theta, b)} for the first likelihood evaluation.
#'   The default (\code{\link{NULL}}) is to accept the internally generated
#'   default.
#'   Use \code{fit = FALSE} to retrieve this default and other optimization
#'   inputs (in particular \code{Y_init}; see Value), which can often be used
#'   to construct a more informative starting point.
#' @param map
#'   A \link{factor} corresponding elementwise to the full parameter vector
#'   \code{c(beta, theta, b)}.
#'   Elements of c(beta, theta, b) corresponding to \code{\link{NA}} in
#'   \code{map} are fixed at their initial values, rather than estimated.
#'   Elements corresponding to a common factor level are constrained to have
#'   a common value during estimation. (However, note that grouping across
#'   segments \code{beta}, \code{theta}, and \code{b} is not implemented.)
#'   \link[=logical]{Logical} values of \code{map} are accepted and coerced
#'   to factor internally with
#'   \code{map <- \link{replace}(\link{gl}(\link{length}(map), 1), map, NA)}.
#' @param append
#'   An expression indicating variables in \code{data_windows}
#'   to be preserved in the returned object for use by methods.
#'   Usage requires that \code{data_windows} is a
#'   \link[=data.frame]{data frame}.
#'   The default (\code{\link{NULL}}) is to preserve nothing.
#'   A dot \samp{.} is to preserve all variables not occurring
#'   in \code{formula_parameters}.
#'   Outside of these two special cases, expressions are evaluated
#'   similarly to argument \code{select} of function \code{\link{subset}}.
#' @param ...
#'   Arguments passed to methods by the generic function.
#'
#' @details
#' Let
#' \code{formula = cbind(time, x) ~ ts}
#' and
#' \code{formula_windows = cbind(start, end) ~ ts_windows}.
#' It is expected that \code{time}, \code{start}, and \code{end}
#' (after possible coercion to \link{numeric}) measure time
#' on the same scale. To be precise, numeric times should have
#' a common unit of measure and, at least within time series,
#' represent displacements from a common reference time.
#' These conditions will always hold if \code{time}, \code{start},
#' and \code{end} all evaluate to \link{Date} or \link{POSIXt} vectors.
#'
#' Plot methods asked to display time on a Date axis rather than
#' a numeric axis will interpret numeric times as numbers of days
#' since \code{1970-01-01 00:00:00}. This interpretation will
#' always be correct if \code{time}, \code{start}, and \code{end}
#' all evaluate to \link{Date} or \link{POSIXt} vectors.
#'
#' When day of week effects are estimated (model$day_of_week > 0),
#' numeric times \emph{must} be interpretable as numbers of days
#' since \code{1970-01-01 00:00:00}, so that time points can be
#' mapped unambiguously to days of week. Furthermore, in this case,
#' \link{time} (after possible coercion to \link{numeric}) is required
#' to be integer-valued with one day spacing in all time series.
#' This means that
#' \code{\link{all.equal}(time, \link{round}(time))}
#' and
#' \code{\link{all}(\link{diff}(\link{round}(time)) == 1)}
#' must return \code{TRUE} in each level of \code{ts}.
#' These conditions ensure that intervals between adjacent
#' time points capture a single day of week.
#'
#' @return
#' If \code{fit = TRUE}, then a \link{list} inheriting from \link{class}
#' \code{"egf"}, with elements:
#' \item{model}{
#'   A copy of the so-named argument.
#' }
#' \item{frame}{
#'   A \link[=data.frame]{data frame} constructed from \code{formula},
#'   \code{data}, and \code{subset}, with variables \code{ts}, \code{time},
#'   and \code{x} specifying time series relevant to the fitted model
#'   (i.e., those with at least one valid fitting window). Variable
#'   \code{window} groups the elements of \code{x} by fitting window.
#'   Rows are ordered by time series and chronologically within time series.
#'   \code{formula} is retained as an \link[=attributes]{attribute}.
#' }
#' \item{frame_windows}{
#'   A \link[=data.frame]{data frame} constructed from \code{formula_windows},
#'   \code{data_windows}, and \code{subset_windows}, with variables \code{ts},
#'   \code{window}, \code{start}, and \code{end} listing fitting windows
#'   for each time series. These intervals may not match the ones supplied
#'   in the function call, which are contracted internally to the narrowest
#'   intervals containing the same set of observations.
#'   Rows are ordered by time series and chronologically within time series.
#'   \code{formula_windows} is retained as an \link[=attributes]{attribute}.
#' }
#' \item{frame_parameters}{
#'   A \link{list} of mixed effects \link[=model.frame]{model frames}
#'   (one for each top level nonlinear model parameter) constructed from
#'   \code{formula_parameters}, \code{data_windows}, and \code{subset_windows}.
#'   \code{frame_parameters[[name]]} retains
#'   \code{\link{terms}(formula_parameters[[name]])} as an
#'   \link[=attributes]{attribute}.
#'   Model frames correspond rowwise to \code{frame_windows}.
#' }
#' \item{frame_append}{
#'   A \link[=data.frame]{data frame} preserving variables
#'   from \code{data_windows} indicated by \code{append}.
#'   This corresponds rowwise to \code{frame_windows}.
#' }
#' \item{priors_top, priors_bottom}{
#'   Named \link{list}s of \code{"\link{egf_prior}"} objects obtained after
#'   processing of \code{formula_priors_top} and \code{formula_priors_bottom}.
#'   \code{priors_top} has one element for each top level nonlinear model
#'   parameter.
#'   \code{priors_bottom} corresponds elementwise to \code{c(beta, theta)}.
#' }
#' \item{control}{
#'   A copy of the so-named argument.
#' }
#' \item{tmb_args}{
#'   A \link{list} of arguments to \code{\link[TMB]{MakeADFun}},
#'   enabling reinitialization of \code{tmb_out}.
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
#'   \code{c(beta, theta, b)} of the first and best likelihood evaluations.
#' }
#' \item{nonrandom}{
#'   An \link{integer} vector indexing segments \code{beta} and \code{theta}
#'   of the full parameter vector \code{c(beta, theta, b)}. These are the
#'   elements that are \emph{not} random effects.
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
#' If \code{fit = FALSE}, then a \link{list} inheriting from
#' \link{class} \code{"egf_no_fit"} containing all of the above
#' except \code{optimizer_out}, \code{best}, \code{sdreport},
#' and \code{nll}. It contains an additional element
#' \code{Y_init}, a \link[=double]{numeric} \link{matrix}
#' supplying naive estimates of each top level nonlinear model
#' parameter, corresponding rowwise to \code{frame_windows}.
#'
#' @export
#' @useDynLib epigrowthfit
egf <- function(model, ...) {
  UseMethod("egf", model)
}

#' @rdname egf
#' @export
#' @importFrom TMB MakeADFun openmp
egf.egf_model <- function(model,
                          formula,
                          formula_windows,
                          formula_parameters = list(),
                          formula_priors_top = list(),
                          formula_priors_bottom = list(),
                          data,
                          data_windows,
                          subset = NULL,
                          subset_windows = NULL,
                          na_action = c("fail", "pass"),
                          na_action_windows = c("fail", "omit"),
                          control = egf_control(),
                          fit = TRUE,
                          se = FALSE,
                          init = NULL,
                          map = NULL,
                          append = NULL,
                          ...) {
  stopifnot(
    inherits(formula, "formula"),
    inherits(formula_windows, "formula")
  )
  if (inherits(formula_parameters, "formula")) {
    stopifnot(length(formula_parameters) == 2L)
  } else {
    stopifnot(
      is.list(formula_parameters),
      vapply(formula_parameters, inherits, FALSE, "formula"),
      lengths(formula_parameters) == 3L
    )
  }
  stopifnot(
    is.list(formula_priors_top),
    vapply(formula_priors_top, inherits, FALSE, "formula"),
    lengths(formula_priors_top) == 3L
  )
  stopifnot(
    is.list(formula_priors_bottom),
    vapply(formula_priors_bottom, inherits, FALSE, "formula"),
    lengths(formula_priors_bottom) == 3L
  )
  if (missing(data)) {
    data <- environment(formula)
  } else {
    stopifnot(is.list(data) || is.environment(data))
  }
  if (missing(data_windows)) {
    data_windows <- environment(formula_windows)
  } else {
    stopifnot(is.list(data_windows) || is.environment(data_windows))
  }
  subset <- substitute(subset)
  subset_windows <- substitute(subset_windows)
  na_action <- match.arg(na_action)
  na_action_windows <- match.arg(na_action_windows)
  stopifnot(inherits(control, "egf_control"))
  stop_if_not_true_false(fit)
  stop_if_not_true_false(se)
  stopifnot(is.numeric(init) || is.null(init))
  if (is.logical(map)) {
    map <- replace(gl(length(map), 1L), map, NA)
  } else {
    stopifnot(is.factor(map) || is.null(map))
  }
  append <- substitute(append)

  formula <- egf_sanitize_formula(formula)
  formula_windows <- egf_sanitize_formula(formula_windows)
  formula_parameters <- egf_sanitize_formula_parameters(formula_parameters, model = model, ignore_intercept = !is.null(init))
  frames <- egf_make_frames(
    model = model,
    formula = formula,
    formula_windows = formula_windows,
    formula_parameters = formula_parameters,
    data = data,
    data_windows = data_windows,
    subset = subset,
    subset_windows = subset_windows,
    na_action = na_action,
    na_action_windows = na_action_windows,
    append = append
  )
  tmb_args <- egf_make_tmb_args(
    model = model,
    frame = frames$frame,
    frame_parameters = frames$frame_parameters,
    control = control,
    fit = fit,
    init = init,
    map = map
  )
  priors_top <- egf_make_priors_top(
    formula_priors_top = formula_priors_top,
    model = model
  )
  priors_bottom <- egf_make_priors_bottom(
    formula_priors_bottom = formula_priors_bottom,
    beta = tmb_args$parameters$beta,
    theta = if (egf_has_random(tmb_args$data)) tmb_args$parameters$theta else numeric(0L)
  )
  tmb_args <- egf_update_tmb_args(
    tmb_args = tmb_args,
    priors_top = priors_top,
    priors_bottom = priors_bottom
  )

  on <- openmp(n = NULL)
  if (on > 0L) {
    openmp(n = control$omp_num_threads)
    on.exit(openmp(n = on))
  }
  tmb_out <- do.call(MakeADFun, tmb_args)
  tmb_out$fn <- egf_patch_fn(tmb_out$fn, inner_optimizer = control$inner_optimizer)
  tmb_out$gr <- egf_patch_gr(tmb_out$gr, inner_optimizer = control$inner_optimizer)
  init <- enum_dupl_names(tmb_out$env$par)
  nonrandom <- seq_along(init)
  if (!is.null(tmb_out$env$random)) {
    nonrandom <- nonrandom[-tmb_out$env$random]
  }

  if (!fit) {
    res <- frames[c("frame", "frame_windows", "frame_parameters")]
    res <- c(res, list(
      model = model,
      tmb_args = tmb_args,
      tmb_out = tmb_out,
      init = init,
      nonrandom = nonrandom,
      Y_init = attr(tmb_args$parameters, "Y_init"),
      call = match.call()
    ))
    class(res) <- c("egf_no_fit", "list")
    return(res)
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

  res <- c(frames, list(
    model = model,
    priors_top = priors_top,
    priors_bottom = priors_bottom,
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
  class(res) <- c("egf", "list")
  res
}
