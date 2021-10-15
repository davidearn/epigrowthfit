#' Fit nonlinear mixed effects models of epidemic growth
#'
#' @description
#' Fits nonlinear mixed effects models of epidemic growth
#' to collections of one or more disease incidence time series.
#'
#' @param model
#'   An \code{"\link{egf_model}"} object specifying a top level
#'   nonlinear model to be estimated.
#' @param formula_ts
#'   A \link{formula} of the form \code{cbind(time, x) ~ ts}
#'   specifying one or more disease incidence time series in long format.
#'   \code{ts} must evaluate to a \link{factor}
#'   (insofar as \code{\link{as.factor}(ts)} is a factor)
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
#'   \code{formula_ts = cbind(time, x) ~ 1} can be supplied
#'   when there is only one time series; it is equivalent
#'   to \code{formula_ts = cbind(time, x) ~ ts} with \code{ts}
#'   evaluating to \code{\link{rep}(factor(1), length(x))}.
#' @param formula_windows
#'   A \link{formula} of the form \code{cbind(start, end) ~ ts}
#'   specifying disjoint fitting windows \code{(start, end]} in long format.
#'   If
#'   \code{formula_ts = cbind(time, x) ~ ts1}
#'   and
#'   \code{formula_windows = cbind(start, end) ~ ts2},
#'   then observation \code{x[i]}
#'   is associated with window \code{(start[j], end[j]]} if and only if
#'   \code{time[i-1] >= start[j]},
#'   \code{time[i] <= end[j]}, and
#'   \code{ts1[i] == ts2[j]}.
#' @param formula_parameters
#'   A \link{list} of \link{formula}e of the form \code{parameter ~ terms}
#'   specifying mixed effects models for top level nonlinear model parameters,
#'   using \pkg{lme4}-like syntax (see \code{?lme4::lmer}).
#'   Alternatively, a formula of the form \code{~terms} to be recycled for
#'   all parameters.
#'   A list of parameters for which formulae may be specified can be retrieved
#'   with \code{\link{egf_get_names_top}}.
#'   Specifically, \code{\link{deparse}(parameter)} must be an element of
#'   \code{\link{egf_get_names_top}(model, link = TRUE)}.
#'   The default for parameters not assigned a formula is \code{~1}.
#' @param formula_priors
#'   A \link{list} of \link{formula}e of the form \code{parameter ~ prior}
#'   defining priors on:\cr
#'   (i) top level nonlinear model parameters,\cr
#'   (ii) fixed effects coefficients and random effect covariance parameters
#'   (elements of segments \code{beta} and \code{theta} of the bottom level
#'   parameter vector \code{c(beta, theta, b)}), or\cr
#'   (iii) random effect covariance matrices
#'   (elements of a \link{list} \code{Sigma} containing the matrices).\cr
#'   \code{prior} must be a \link{call} to a \link[=egf_prior]{prior} function
#'   with arguments specifying suitable hyperparameters.
#'   In case (i),
#'   \code{deparse(parameter)} must be an element of
#'   \code{\link{egf_get_names_top}(model, link = TRUE)},
#'   and hyperparameters supplied on the right hand side must have length 1.
#'   In cases (ii) and (iii),
#'   \code{parameter} must be \code{beta}, \code{theta}, or \code{Sigma}
#'   or a call to \code{\link{[}} or \code{\link{[[}}
#'   subsetting \code{beta}, \code{theta}, or \code{Sigma}
#'   (e.g., \code{beta[index]}, where \code{index} is a valid index vector
#'   for \code{beta}),
#'   and hyperparameters are recycled to the length of the indicated subset.
#'   All expressions \code{prior} and \code{index} are evaluated in the
#'   corresponding formula environment.
#' @param data_ts,data_windows
#'   \link[=data.frame]{Data frame}s, \link{list}s, or \link{environment}s
#'   to be searched for variables named in the corresponding formulae and
#'   subset expressions. (\code{data_windows} is also searched for variables
#'   named in \code{formula_parameters}.) Formula environments are searched
#'   for variables not found here.
#' @param subset_ts,subset_windows
#'   Expressions to be evaluated in the corresponding data frame.
#'   The result should be a valid index vector for the rows of the data frame
#'   (see \code{\link{[.data.frame}}).
#'   Rows that are not indexed are discarded.
#'   Rows that are indexed are filtered further
#'   (e.g., time series with zero associated fitting windows are discarded
#'   regardless of \code{subset}).
#'   The default (\code{\link{NULL}}) is to preserve all rows for further
#'   filtering.
#' @param na_action_ts
#'   A \link{character} string affecting the handling of \code{\link{NA}}
#'   in \code{x} if \code{formula_ts = cbind(time, x) ~ ts}.
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
#'   A \link{logical} flag. If \code{FALSE}, then \code{egf} returns early
#'   (\emph{before} fitting) with a partial model object.
#' @param se
#'   A \link{logical} flag.
#'   If \code{TRUE}, then the Hessian matrix of the negative log likelihood
#'   function is computed and inverted to approximate the joint covariance
#'   matrix of segments \code{beta} and \code{theta} of the bottom level
#'   parameter vector \code{c(beta, theta, b)}.
#'   In addition, the standard errors of the fitted values of all top level
#'   nonlinear model parameters are computed approximately using the delta
#'   method.
#'   Computations are preserved in the model object for reuse by methods.
#' @param init
#'   A \link{numeric} vector to be used as the bottom level parameter vector
#'   \code{c(beta, theta, b)} for the first likelihood evaluation.
#'   The default (\code{\link{NULL}}) is to accept the internally generated
#'   default, which is a zero vector except for the \code{"(Intercept)"}
#'   coefficients in \code{beta}.
#' @param map
#'   A \link{factor} corresponding elementwise to the bottom level parameter
#'   vector \code{c(beta, theta, b)}.
#'   Elements of \code{c(beta, theta, b)} corresponding to \code{\link{NA}}
#'   in \code{map} are fixed at their initial values, rather than estimated.
#'   Elements corresponding to a common factor level (within segments) are
#'   constrained to have a common value during estimation.
#'   Alternatively, \code{map} can be an index vector for the elements of
#'   \code{c(beta, theta, b)}.
#'   In this case, the indexed elements are fixed at their initial values.
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
#' If
#' \code{formula_ts = cbind(time, x) ~ ts1}
#' and
#' \code{formula_windows = cbind(start, end) ~ ts2},
#' then
#' it is expected that \code{time}, \code{start}, and \code{end}
#' (after possible coercion to \link{numeric}) measure time
#' on the same scale. To be precise, numeric times should have
#' a common unit of measure and, at least within time series,
#' represent displacements from a common reference time.
#' These conditions will always hold if \code{time}, \code{start},
#' and \code{end} all evaluate to \link{Date} or \link{POSIXt}
#' vectors.
#'
#' When day of week effects are estimated (model$day_of_week > 0),
#' numeric times must be interpretable as numbers of days since
#' \code{1970-01-01 00:00:00}, so that time points can be mapped
#' unambiguously to days of week. Furthermore, in this case,
#' \code{time} (after possible coercion to \link{numeric})
#' is required to be integer-valued with one day spacing in all
#' time series. This means that
#' \code{\link{all.equal}(time, \link{round}(time))}
#' and
#' \code{\link{all}(\link{diff}(\link{round}(time)) == 1)}
#' must both be \code{TRUE} in each level of \code{ts}.
#' These conditions ensure that intervals between successive
#' time points each represent exactly one day of week.
#'
#' @return
#' A \link{list} inheriting from \link{class} \code{"egf"}
#' (or \code{"egf_no_fit"} if \code{fit = FALSE}).
#' See topic \code{\link{egf-class}} for a description of the list elements.
#' Only developers should need to access the list directly;
#' typical users can rely on methods to retrieve information
#' about the estimated (or to-be-estimated) model.
#' Available methods can be queried by running, e.g.,
#' \code{\link{methods}(class = "egf")}.
#' Links to method help pages can be obtained by running, e.g.,
#' \code{??`\\.egf$`} in an interactive \R session.
#'
#' @examples
#' ## Simulate 'N' incidence time series exhibiting exponential growth
#' set.seed(180149L)
#' N <- 10L
#' time <- seq.int(0, 40, 1)
#' r <- rlnorm(N, -3.2, 0.2)
#' c0 <- rlnorm(N, 6, 0.2)
#' f <- function(time, r, c0) {
#'   lambda <- diff(exp(log(c0) + r * time))
#'   c(NA, rpois(lambda, lambda))
#' }
#' data_ts <- data.frame(
#'   country = gl(N, length(time), labels = LETTERS[seq_len(N)]),
#'   time = rep.int(time, N),
#'   x = unlist(Map(f, time = list(time), r = r, c0 = c0))
#' )
#'
#' ## Define fitting windows (here, two per time series)
#' data_windows <- data.frame(
#'   country = gl(N, 1L, 2L * N, labels = LETTERS[seq_len(N)]),
#'   wave = gl(2L, 10L),
#'   start = c(
#'     sample(seq.int(0, 5, 1), size = N, replace = TRUE),
#'     sample(seq.int(20, 25, 1), size = N, replace = TRUE)
#'   ),
#'   end = c(
#'     sample(seq.int(15, 20, 1), size = N, replace = TRUE),
#'     sample(seq.int(35, 40, 1), size = N, replace = TRUE)
#'   )
#' )
#'
#' ## Estimate the generative model
#' root <- system.file(package = "epigrowthfit", mustWork = TRUE)
#' path_to_cache <- file.path(root, "exdata", "egf.rds")
#' if (file.exists(path_to_cache)) {
#'   object <- readRDS(path_to_cache)
#' } else {
#'   object <- egf(
#'     model = egf_model(curve = "exponential", family = "pois"),
#'     formula_ts = cbind(time, x) ~ country,
#'     formula_windows = cbind(start, end) ~ country,
#'     formula_parameters = ~(1 | country:wave),
#'     formula_priors = list(Sigma ~ LKJ(eta = 2)),
#'     data_ts = data_ts,
#'     data_windows = data_windows,
#'     se = TRUE
#'   )
#'   dir.create(dirname(path_to_cache), showWarnings = FALSE)
#'   saveRDS(object, file = path_to_cache)
#' }
#'
#' @export
#' @useDynLib epigrowthfit
egf <- function(model, ...) {
  UseMethod("egf", model)
}

#' @rdname egf
#' @export
#' @importFrom TMB MakeADFun openmp sdreport
egf.egf_model <- function(model,
                          formula_ts,
                          formula_windows,
                          formula_parameters = list(),
                          formula_priors = list(),
                          data_ts,
                          data_windows,
                          subset_ts = NULL,
                          subset_windows = NULL,
                          na_action_ts = c("fail", "pass"),
                          na_action_windows = c("fail", "omit"),
                          control = egf_control(),
                          fit = TRUE,
                          se = FALSE,
                          init = NULL,
                          map = NULL,
                          append = NULL,
                          ...) {
  stopifnot(
    inherits(formula_ts, "formula"),
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
    is.list(formula_priors),
    vapply(formula_priors, inherits, FALSE, "formula"),
    lengths(formula_priors) == 3L
  )
  if (missing(data_ts)) {
    data_ts <- environment(formula_ts)
  } else {
    stopifnot(is.list(data_ts) || is.environment(data_ts))
  }
  if (missing(data_windows)) {
    data_windows <- environment(formula_windows)
  } else {
    stopifnot(is.list(data_windows) || is.environment(data_windows))
  }
  subset_ts <- substitute(subset_ts)
  subset_windows <- substitute(subset_windows)
  na_action_ts <- match.arg(na_action_ts)
  na_action_windows <- match.arg(na_action_windows)
  stopifnot(inherits(control, "egf_control"))
  stop_if_not_true_false(fit)
  stop_if_not_true_false(se)
  stopifnot(is.numeric(init) || is.null(init))
  append <- substitute(append)

  names_parameters <- egf_get_names_top(model, link = TRUE)

  formula_ts <- egf_sanitize_formula(formula_ts)
  formula_windows <- egf_sanitize_formula(formula_windows)
  formula_parameters <- egf_sanitize_formula_parameters(
    formula_parameters = formula_parameters,
    names_parameters = names_parameters,
    check_intercept = is.null(init)
  )
  frame <- egf_make_frame(
    model = model,
    formula_ts = formula_ts,
    formula_windows = formula_windows,
    formula_parameters = formula_parameters,
    data_ts = data_ts,
    data_windows = data_windows,
    subset_ts = subset_ts,
    subset_windows = subset_windows,
    na_action_ts = na_action_ts,
    na_action_windows = na_action_windows,
    append = append
  )

  env <- new.env(parent = emptyenv())
  tmb_args <- egf_tmb_make_args(
    model = model,
    frame = frame,
    control = control,
    init = init,
    map = map,
    env = env
  )

  priors <- egf_make_priors(
    formula_priors = formula_priors,
    top = list(
      names = names_parameters,
      family = "norm"
    ),
    beta = list(
      length = env$len[["beta"]],
      family = "norm"
    ),
    theta = list(
      length = env$len[["theta"]],
      family = "norm"
    ),
    Sigma = list(
      length = length(tmb_args$data$block_rows),
      family = c("lkj", "wishart", "invwishart"),
      rows = tmb_args$data$block_rows
    )
  )
  tmb_args$data <- egf_tmb_update_data(tmb_args$data, priors = priors)

  tmb_out <- do.call(MakeADFun, tmb_args)
  tmb_out$fn <- egf_patch_fn(tmb_out$fn, inner_optimizer = control$inner_optimizer)
  tmb_out$gr <- egf_patch_gr(tmb_out$gr, inner_optimizer = control$inner_optimizer)

  res <- list(
    model = model,
    frame = frame,
    priors = priors,
    control = control,
    tmb_out = tmb_out,
    optimizer_out = NULL,
    init = tmb_out$env$par,
    best = NULL,
    random = tmb_out$env$lrandom(),
    value = NULL,
    gradient = NULL,
    hessian = NULL,
    sdreport = NULL,
    effects = env$effects,
    contrasts = env$contrasts,
    Y0 = env$Y0,
    call = match.call()
  )

  if (!fit) {
    class(res) <- "egf_no_fit"
    return(res)
  }

  on <- openmp(n = NULL)
  if (on > 0L) {
    openmp(n = control$omp_num_threads)
    on.exit(openmp(n = on))
  }
  optimizer <- control$optimizer$f
  optimizer_args <- c(
    tmb_out[c("par", "fn", "gr")],
    control$optimizer["control"],
    control$optimizer[["args"]]
  )
  res$optimizer_out <- do.call(optimizer, optimizer_args)

  res$best <- tmb_out$env$last.par.best
  res$value <- as.numeric(tmb_out$env$value.best)
  if (se) {
    res$sdreport <- try(sdreport(tmb_out, par.fixed = res$best[!res$random], getReportCovariance = FALSE))
  }
  if (inherits(res$sdreport, "sdreport")) {
    res$gradient <- res$sdreport$gradient.fixed
    res$hessian <- res$sdreport$pdHess
  } else {
    res$gradient <- tmb_out$gr(res$best[!res$random])
    res$hessian <- NA
  }

  class(res) <- "egf"
  res
}

#' Fitted and unfitted egf objects
#'
#' An object returned by \code{\link{egf}},
#' inheriting from class \code{"egf"} if the nonlinear mixed effects model
#' that it specifies has been fitted and class \code{"egf_no_fit"} otherwise.
#'
#' @return
#' Legitimate \code{"egf"} and \code{"egf_no_fit"} objects are lists
#' with elements:
#' \item{model}{
#'   A copy of the so-named argument.
#' }
#' \item{frame}{
#'   A recursive list of data frames.
#' }
#' \item{priors}{
#'   A recursive list of \code{"\link{egf_prior}"} objects,
#'   of the form \code{list(top, bottom = list(beta, theta, Sigma))}.
#' }
#' \item{control}{
#'   A copy of the so-named argument.
#' }
#' \item{tmb_out}{
#'   The list output of \code{\link[TMB]{MakeADFun}}.
#' }
#' \item{optimizer_out}{
#'   The \link{list} output of the optimizer specified by
#'   \code{control$optimizer}.
#' }
#' \item{init, best}{
#'
#' }
#' \item{random}{
#'
#' }
#' \item{value, gradient}{
#'   Numeric vectors giving the value and gradient of the negative
#'   log (marginal) likelihood function.
#' }
#' \item{hessian}{
#'   A logical flag indicating whether the Hessian matrix of the negative
#'   log (marginal) likelihood function, evaluated at \code{best[!random]},
#'   is positive definite.
#'   \code{TRUE} means the matrix is positive definite.
#'   \code{NA} means the matrix has not been computed,
#'   either because \code{se = FALSE} or because an error
#'   was thrown during computation by \code{\link{sdreport}}.
#' }
#' \item{sdreport}{
#'   If \code{se = TRUE}, then a list inheriting from class \code{"sdreport"},
#'   resulting from \code{\link{sdreport}(tmb_out)}. Otherwise, \code{NULL}.
#' }
#' \item{effects}{
#'
#' }
#' \item{contrasts}{
#'
#' }
#' \item{Y0}{
#'
#' }
#' \item{call}{
#'   The \link{call} to \code{\link{egf}},
#'   enabling updates to the object via \code{\link{update}}.
#' }
#' \code{optimizer_out}, \code{best},
#' \code{value}, \code{gradient}, \code{hessian}, and \code{sdreport}
#' are \code{NULL} in \code{"egf_no_fit"} objects.
#'
#' @name egf-class
NULL
