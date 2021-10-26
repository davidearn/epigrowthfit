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
#'   \code{\link{egf_get_names_top}(model)}.
#'   The default for parameters not assigned a formula is \code{~1}.
#' @param formula_priors
#'   A \link{list} of \link{formula}e of the form \code{parameter ~ prior}
#'   defining priors on:\cr
#'   (i) top level nonlinear model parameters,\cr
#'   (ii) fixed effects coefficients and random effect covariance parameters
#'   (elements of segments \code{beta} and \code{theta} of the bottom level
#'   parameter vector, or\cr
#'   (iii) random effect covariance matrices
#'   (elements of a \link{list} \code{Sigma} containing the matrices).\cr
#'   \code{prior} must be a \link{call} to a \link[=egf_prior]{prior} function
#'   with arguments specifying suitable hyperparameters.
#'   In case (i),
#'   \code{deparse(parameter)} must be an element of
#'   \code{\link{egf_get_names_top}(model)},
#'   and hyperparameters supplied on the right hand side must have length 1.
#'   In cases (ii) and (iii),
#'   \code{parameter} must be \code{beta}, \code{theta}, or \code{Sigma}
#'   or a call to \code{\link{[}} or \code{\link{[[}} referring to a subset
#'   or element of \code{beta}, \code{theta}, or \code{Sigma}
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
#'   parameter vector.
#'   In addition, the standard errors of the fitted values of all top level
#'   nonlinear model parameters are computed approximately using the delta
#'   method.
#'   Computations are preserved in the model object for reuse by methods.
#' @param init
#'   A named \link{list} of \link{numeric} vectors with possible elements
#'   \code{beta}, \code{theta}, and \code{b}, specifying values to be used
#'   in the first likelihood evaluation for the so-named segments of the
#'   bottom level parameter vector.
#'   The default value of each segment is a zero vector, with the exception
#'   that \code{"(Intercept)"} coefficients in \code{beta} have default values
#'   computed from supplied time series.
#'   To specify only a subset of a segment, use \code{\link{NA}} to indicate
#'   elements that should retain their default value.
#' @param map
#'   A named list of \link{factor}s with possible elements \code{beta},
#'   \code{theta}, and \code{b}, each as long as the so-named segment
#'   of the bottom level parameter.
#'   Elements of segment \code{name} indexed by \code{is.na(map[["name"]])},
#'   are fixed at their initial values, rather than estimated, and elements
#'   corresponding to a common factor level are are constrained to have a
#'   common value during estimation.
#'   \code{map[["name"]]} can be an index vector for segment \code{name},
#'   instead of a factor. In this case, the elements of segment \code{name}
#'   indexed by \code{map[["name"]]} are fixed at their initial values.
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
#' Users attempting to set arguments \code{formula_priors}, \code{init}, and
#' \code{map} should know the structure of the bottom level parameter vector.
#' It is described under topic \code{\link{egf-class}}.
#'
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
#' must both be \code{TRUE} in each level of \code{ts1}.
#' These conditions ensure that intervals between successive
#' time points each represent exactly one day of week.
#'
#' @return
#' A \link{list} inheriting from \link{class} \code{"egf"}
#' or \code{"egf_no_fit"} depending on the value of \code{fit}.
#' See topic \code{\link{egf-class}} for class documentation.
#'
#' @examples
#' ## Simulate 'N' incidence time series exhibiting exponential growth
#' set.seed(180149L)
#' N <- 10L
#' time <- seq.int(0, 40, 1)
#' mu <- c(-3.2, 6)
#' sigma <- c(0.2, 0.2)
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
#' m1 <- egf_cache("egf-1.rds", {
#'   egf(
#'     model = egf_model(curve = "exponential", family = "pois"),
#'     formula_ts = cbind(time, x) ~ country,
#'     formula_windows = cbind(start, end) ~ country,
#'     formula_parameters = ~(1 | country:wave),
#'     data_ts = data_ts,
#'     data_windows = data_windows,
#'     se = TRUE
#'   )
#' })
#'
#' ## Re-estimate the generative model with:
#' ## * Gaussian prior on 'beta[1L]'
#' ## * LKJ prior on all random effect covariance matrices (though here there is just one)
#' ## * initial value of 'theta' set explicitly
#' ## * 'theta[3L]' fixed at initial value
#' m2 <- egf_cache("egf-2.rds", {
#'   update(m1,
#'     formula_priors = list(
#'       beta[1L] ~ Normal(mu = -3, sigma = 1),
#'       Sigma ~ LKJ(eta = 2)
#'     ),
#'     init = list(theta = c(log(0.5), log(0.5), 0)),
#'     map = list(theta = 3L)
#'   )
#' })
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
                          init = list(),
                          map = list(),
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
  stopifnot(
    inherits(control, "egf_control"),
    is_true_or_false(fit),
    is_true_or_false(se),
    is.list(init),
    is.list(map)
  )
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
  res$value <- as.double(tmb_out$env$value.best)
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

#' Object of class "egf" or "egf_no_fit"
#'
#' An object returned by \code{\link{egf}}, inheriting from class \code{"egf"}
#' or \code{"egf_no_fit"} depending on whether the nonlinear mixed effects model
#' specified in the call was actually fit.
#'
#' @details
#' Only developers should need to access the list directly,
#' e.g., using \code{$}.
#' Typical users can rely on methods to get information about the estimated
#' (or to-be-estimated) model.
#' See Examples for ways to query available methods and their documentation.
#'
#' The estimated (or to-be-estimated) model is specified by a bottom level
#' parameter vector that is the concatenation of three segments:
#' \describe{
#' \item{beta}{
#'   The result of \code{unlist(lbeta)}, where \code{lbeta} is a list of
#'   numeric vectors of fixed effects coefficients, with one vector for
#'   each top level nonlinear model parameter.
#'   The order of top level parameters is given by
#'   \code{\link{egf_get_names_top}(model)}.
#' }
#' \item{theta}{
#'   The result of \code{unlist(ltheta)}, where \code{ltheta} is a list of
#'   numeric vectors of random effect covariance parameters, with one vector
#'   for each distinct random effects term in \code{formula_parameters}.
#'   Each vector parametrizes a random effect covariance matrix via
#'   \code{\link{theta2cov}} and its inverse \code{\link{cov2theta}}.
#'   (The list \code{Sigma} described under \code{\link{egf}} argument
#'   \code{formula_priors} is precisely
#'   \code{\link{lapply}(ltheta, \link{theta2cov})}.)
#' }
#' \item{b}{
#'   The result of \code{unlist(lb)}, where \code{lb} is a list of numeric
#'   matrices of scaled random effects coefficients, corresponding elementwise
#'   to \code{ltheta}.
#'   The columns of \code{lb[[i]]} (one per level of the grouping variable)
#'   are interpreted as samples from a zero mean, unit variance multivariate
#'   normal distribution with correlation matrix
#'   \code{\link{cov2cor}(\link{theta2cov}(ltheta[[i]]))}.
#' }
#' }
#' When elements of this vector are mapped via \code{\link{egf}} argument
#' \code{map}, likelihood is defined as a function of the condensed vector
#' that excludes mapped elements.
#'
#' A number of methods are available to allow users to investigate the
#' structure of each of the three segments, with or without fitting a model;
#' see
#' \code{\link[=coef.egf]{coef}},
#' \code{\link[=fixef.egf]{fixef}}, and
#' \code{\link[=ranef.egf]{ranef}}.
#'
#' @return
#' A legitimate \code{"egf"} or \code{"egf_no_fit"} object is a list
#' with elements:
#' \item{model}{
#'   A copy of the so-named argument of \code{\link{egf}}.
#' }
#' \item{frame}{
#'   A list of the form \code{list(ts, windows, parameters, append)}.
#'   \code{ts} and \code{windows} are data frames preserving time series
#'   and fitting window endpoints.
#'   \code{parameters} is a list of mixed effects model frames,
#'   with one element for each top level nonlinear model parameter.
#'   \code{append} is a data frame preserving additional variables
#'   specified in \code{call$append}.
#'   \code{windows}, the model frames listed in \code{parameters},
#'   and \code{append} all correspond rowwise.
#' }
#' \item{priors}{
#'   A list of the form \code{list(top, bottom = list(beta, theta, Sigma))},
#'   where \code{top}, \code{beta}, \code{theta}, and \code{Sigma} are all
#'   lists of \code{"\link{egf_prior}"} objects.
#' }
#' \item{control}{
#'   A copy of the so-named argument of \code{\link{egf}}.
#' }
#' \item{tmb_out}{
#'   The list output of \code{\link[TMB]{MakeADFun}}.
#'   This contains an \link{environment} \code{env} whose objects
#'   are updated with each evaluation of the objective function.
#' }
#' \item{optimizer_out}{
#'   The list output of the optimizer specified by \code{control$optimizer}.
#' }
#' \item{init, best}{
#'   Numeric vectors giving the values of the condensed bottom level parameter
#'   vector used in the first and maximal likelihood evaluations.
#'   The \code{\link{names}} attribute of each vector groups the elements by
#'   segment; see Details.
#' }
#' \item{random}{
#'   A logical vector indexing the elements of the condensed bottom level
#'   parameter vector that are \emph{not} arguments of the negative log
#'   \emph{marginal} likelihood function.
#'   It indexes all elements of segment \code{b} (random effects) and,
#'   if \code{control$profile = TRUE}, all elements of segment \code{beta};
#'   see Details.
#' }
#' \item{value, gradient}{
#'   Numeric vectors giving the value and gradient of the negative log
#'   marginal likelihood function at \code{best[!random]}.
#'   \code{value} can be extracted using \code{\link[=logLik.egf]{logLik}}.
#' }
#' \item{hessian}{
#'   A logical flag indicating whether the Hessian matrix of the negative log
#'   (marginal) likelihood function is positive definite at \code{best[!random]}.
#'   \code{NA} means the matrix has not been computed
#'   (\code{sdreport} is not an \code{"sdreport"} object).
#' }
#' \item{sdreport}{
#'   If \code{call} contains \code{se = TRUE},
#'   then the result of \code{\link{try}(\link{sdreport}(tmb_out))}.
#'   Otherwise, \code{NULL}.
#' }
#' \item{effects}{
#'   A list of the form \code{list(beta, b)} with \code{beta} and \code{b}
#'   data frames preserving interpretive information about the so-named
#'   segments of the bottom level parameter vector.
#' }
#' \item{contrasts}{
#'   A list of the form \code{list(X, Z)}, with \code{X} and \code{Z}
#'   lists preserving the contrasts used to construct the fixed and
#'   random effects \link[=model.matrix]{design} matrices.
#' }
#' \item{call}{
#'   The \link{call} to \code{\link{egf}},
#'   enabling updates to the object via \code{\link{update}}.
#' }
#' \code{optimizer_out}, \code{best}, \code{value}, \code{gradient},
#' \code{hessian}, and \code{sdreport} have the value \code{NULL} in
#' objects of class \code{"egf_no_fit"}.
#'
#' @examples
#' methods(class = "egf")
#' help.search("\\.egf$", fields = "alias", package = "epigrowthfit")
#' ## less verbosely: alias??`\\\\.egf$`
#'
#' methods(class = "egf_no_fit")
#' help.search("\\.egf_no_fit$", fields = "alias", package = "epigrowthfit")
#' ## less verbosely: alias??`\\\\.egf_no_fit$`
#'
#' @name egf-class
NULL
