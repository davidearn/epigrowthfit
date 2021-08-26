#' Define an epidemic model
#'
#' Sets flags defining the top level nonlinear model of epidemic growth
#' to be estimated by \code{\link{egf}}.
#'
#' @param curve
#'   A \link{character} string specifying a model for expected cumulative
#'   disease incidence as a function of time.
#' @param excess
#'   A \link{logical} flag. If \code{TRUE}, then a constant baseline
#'   mortality rate is estimated. Set to \code{TRUE} if what is
#'   observed is multiple causes mortality rather than disease mortality
#'   or disease incidence.
#' @param family
#'   A \link{character} string specifying a family of discrete probability
#'   distributions assigned to observations of disease incidence.
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
egf_model <- function(curve = c("logistic", "exponential", "subexponential", "gompertz", "richards"),
                      excess = FALSE,
                      family = c("nbinom", "pois"),
                      day_of_week = FALSE) {
  curve <- match.arg(curve)
  stop_if_not_true_false(excess)
  family <- match.arg(family)
  stop_if_not_true_false(day_of_week, allow_numeric = TRUE)
  day_of_week <- as.integer(day_of_week)
  day_of_week <- (day_of_week > 0L) * (1L + (day_of_week - 1L) %% 7L) # coercion to '0:7'

  res <- list(
    curve = curve,
    excess = excess,
    family = family,
    day_of_week = day_of_week
  )
  class(res) <- c("egf_model", "list")
  res
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
#'   satisfied by the nonlinear mixed effects models fit by \code{\link{egf}}.
#' @param sparse_X
#'   A \link{logical} flag. If \code{TRUE}, then the fixed effects
#'   \link[=model.matrix]{design matrix} is constructed in
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
#' \item{2}{
#'   All negative log likelihood terms are printed.
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
  stopifnot(inherits(optimizer, "egf_optimizer"))
  if (inherits(inner_optimizer, "egf_inner_optimizer")) {
    inner_optimizer <- list(inner_optimizer)
  } else {
    stopifnot(
      is.list(inner_optimizer),
      length(inner_optimizer) > 0L,
      vapply(inner_optimizer, inherits, FALSE, "egf_inner_optimizer")
    )
  }
  stop_if_not_true_false(trace, allow_numeric = TRUE)
  trace <- min(2L, max(0L, as.integer(trace))) # coercion to '0:2'
  stop_if_not_true_false(profile)
  stop_if_not_true_false(sparse_X)
  stop_if_not_integer(omp_num_threads, "positive")
  omp_num_threads <- as.integer(omp_num_threads)

  res <- list(
    optimizer = optimizer,
    inner_optimizer = inner_optimizer,
    trace = trace,
    profile = profile,
    sparse_X = sparse_X,
    omp_num_threads = omp_num_threads
  )
  class(res) <- c("egf_control", "list")
  res
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
#'   To align the default behaviour of \code{newton} with that of \code{optim},
#'   \code{trace = 0L} is set if not specified in \code{args}.
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
  optimizer <- f
  stopifnot(
    is.list(args),
    is.list(control)
  )
  if (identical(f, optim)) {
    if (is.null(args$method)) {
      args$method <- "BFGS"
      warning("'optim' argument 'method' not specified, using ", dQuote(args$method), ".")
    } else {
      args$method <- match.arg(args$method, eval(formals(optim)$method))
    }
  } else if (identical(f, nlminb)) {
    f <- function(par, fn, gr, control, ...) {
      res <- nlminb(start = par, objective = fn, gradient = gr, control = control, ...)
      res["value"] <- res["objective"]
      res["objective"] <- NULL
      res
    }
  } else if (identical(f, nlm)) {
    f <- function(par, fn, gr, control, ...) {
      res <- nlm(f = structure(fn, gradient = gr), p = par, ...)
      m <- match(c("estimate", "minimum", "code"), names(res), 0L)
      names(res)[m] <- c("par", "value", "convergence")
      res["message"] <- list(NULL)
      res
    }
  } else {
    stopifnot(
      is.function(f),
      length(nf <- names(formals(f))) >= 4L,
      nf[1:3] != "...",
      "control" %in% nf[-(1:3)]
    )
    e <- quote(f(c(1, 1), function(x) sum(x^2), function(x) 2 * x))
    f_out <- tryCatch(eval(e),
      error = function(cond) {
        stop(wrap(
          "Unable to validate 'f' due to following error in test ", sQuote(deparse(e)), ":\n\n",
          conditionMessage(cond)
        ))
      }
    )
    required <- c("par", "value", "convergence", "message")
    if (!(is.list(f_out) && all(required %in% names(f_out)))) {
      stop(wrap(
        "'f' must return a list with elements ",
        paste(sQuote(required), collapse = ", "),
        " but _does not_ for test ", sQuote(deparse(e)), "."
      ))
    }
    f <- function(par, fn, gr, control, ...) {
      optimizer(par, fn, gr, control = control, ...)
    }
  }
  if (!is.null(names(args))) {
    reserved <- c("par", "fn", "gr", "control", "...", names(formals(optimizer))[1:3])
    args <- args[setdiff(names(args), reserved)]
  }

  res <- list(f = f, args = args, control = control)
  class(res) <- c("egf_optimizer", "list")
  res
}

#' @rdname egf_optimizer
#' @export
#' @importFrom stats optim
#' @importFrom TMB newton
egf_inner_optimizer <- function(f = newton, args = list(), control = list()) {
  stopifnot(
    is.list(args),
    is.list(control)
  )
  if (identical(f, newton)) {
    method <- "newton"
    if (!is.null(names(args))) {
      reserved <- c("par", "fn", "gr", "he", "env", "...")
      args <- args[setdiff(names(args), reserved)]
    }
    if (!"trace" %in% names(args)) {
      args[["trace"]] <- 0L
    }
    control <- args
  } else if (identical(f, optim)) {
    if (is.null(args$method)) {
      method <- "BFGS"
      warning("'optim' argument 'method' not specified, using ", dQuote(method), ".")
    } else {
      method <- match.arg(args$method, eval(formals(optim)$method))
    }
  } else {
    stop("'f' is currently restricted to 'TMB::newton' and 'stats::optim'.")
  }

  res <- list(method = method, control = control)
  class(res) <- c("egf_inner_optimizer", "list")
  res
}

#' Define a parallelization method
#'
#' Defines instructions for parallelization by linking a method with options.
#'
#' @param method
#'   A \link{character} string indicating a method of parallelization.
#'   \code{"\link[=lapply]{serial}"} indicates no parallelization.
#'   \code{"\link[parallel:mclapply]{multicore}"} indicates \R level
#'   forking. It is intended for use from a terminal rather than a GUI.
#'   On Windows, \code{"multicore"} is equivalent to \code{"serial"}.
#'   \code{"\link[parallel:parLapply]{snow}"} indicates socket clusters.
#'   \code{"snow"} is supported on both Unix-alikes and Windows.
#' @param outfile
#'   A \link{character} string indicating a file path where console output
#'   should be diverted. An empty string indicates no diversion.
#'   If \code{method = "snow"}, then diversion may be necessary to view output.
#' @param cores
#'   A positive integer indicating a number of threads/processes
#'   to fork/spawn when \code{parallel != "serial"}. The maximum
#'   is typically \code{\link[parallel]{detectCores}(TRUE, FALSE)}.
#' @param args
#'   A \link{list} of optional arguments to
#'   \code{\link[parallel]{mclapply}} (\code{method = "multicore"}) or
#'   \code{\link[parallel]{makePSOCKcluster}} (\code{method = "snow"}).
#' @param cl
#'   An existing \link[parallel:makePSOCKcluster]{socket cluster}
#'   (\code{method = "snow"}).
#'   The default (\code{\link{NULL}}) is to create a new clusters
#'   as necessary and terminate them upon job completion.
#'   (If non-\code{\link{NULL}}, then \code{outfile}, \code{cores},
#'   and \code{options} are ignored.)
#'
#' @details
#' For general information about parallelism in \R,
#' see \code{\link{vignette}("parallel", "parallel")}.
#'
#' @return
#' A \link{list} inheriting from \link{class} \code{"egf_parallel"}
#' containing the arguments (after possible matching and coercion).
#'
#' @export
#' @importFrom TMB openmp
egf_parallel <- function(method = c("serial", "multicore", "snow"),
                         outfile = nullfile(),
                         cores = getOption("egf.cores", 1L),
                         args = list(),
                         cl = NULL) {
  method <- match.arg(method)
  stop_if_not_string(outfile)
  if (method == "serial") {
    cores <- 1L
    args <- list()
    cl <- NULL
  } else if (method == "multicore" || (method == "snow" && is.null(cl))) {
    stop_if_not_integer(cores, "positive")
    cores <- as.integer(cores)
    stopifnot(is.list(args))
    if (method == "multicore") {
      args$mc.cores <- cores
    } else {
      args$names <- cores
      args$outfile <- outfile
    }
  } else {
    stopifnot(inherits(cl, "SOCKcluster"))
    cores <- length(cl)
    args <- list()
  }

  res <- list(
    method = method,
    outfile = outfile,
    cores = cores,
    args = args,
    cl = cl
  )
  class(res) <- c("egf_parallel", "list")
  res
}

#' Get top level nonlinear model parameter names
#'
#' Retrieves the names used internally for top level nonlinear model
#' parameters.
#'
#' @param object
#'   An \R object specifying a top level nonlinear model,
#'   or \code{\link{NULL}}.
#' @param link
#'   A \link{logical} flag. If \code{TRUE},
#'   then \code{"link(name)"} is returned instead of \code{"name"}.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A \link{character} vector giving the subset of names relevant to
#' \code{object}, or the complete set if \code{object = \link{NULL}}.
#'
#' @examples
#' egf_get_names_top(NULL, link = FALSE)
#' egf_get_names_top(NULL, link = TRUE)
#'
#' model <- egf_model()
#' egf_get_names_top(model, link = FALSE)
#' egf_get_names_top(model, link = TRUE)
#'
#' @export
egf_get_names_top <- function(object, ...) {
  UseMethod("egf_get_names_top", object)
}

#' @rdname egf_get_names_top
#' @export
egf_get_names_top.default <- function(object, link = TRUE, ...) {
  stopifnot(is.null(object))
  names_top <- c("r", "alpha", "c0", "tinfl", "K",
                 "p", "a", "b", "disp", paste0("w", 1:6))
  if (link) {
    return(egf_link_add(names_top))
  }
  names_top
}

#' @rdname egf_get_names_top
#' @export
egf_get_names_top.egf_model <- function(object, link = TRUE, ...) {
  names_top <- switch(object$curve,
    exponential    = c("r", "c0"),
    subexponential = c("alpha", "c0", "p"),
    gompertz       = c("alpha", "tinfl", "K"),
    logistic       = c("r", "tinfl", "K"),
    richards       = c("r", "tinfl", "K", "a")
  )
  if (object$excess) {
    names_top <- c(names_top, "b")
  }
  if (object$family == "nbinom") {
    names_top <- c(names_top, "disp")
  }
  if (object$day_of_week > 0L) {
    names_top <- c(names_top, paste0("w", 1:6))
  }
  if (link) {
    return(egf_link_add(names_top))
  }
  names_top
}

#' @rdname egf_get_names_top
#' @export
egf_get_names_top.egf <- function(object, link = TRUE, ...) {
  egf_get_names_top(object$model, link = link)
}

#' @export
egf_get_names_top.tmb_data <- function(object, link = TRUE, ...) {
  names_top <- levels(object$X_info$par)
  if (link) {
    return(names_top)
  }
  egf_link_remove(names_top)
}

#' Check for random effects
#'
#' Determines whether an object specifies a random effects model.
#'
#' @param object An \R object.
#'
#' @return
#' \code{TRUE} or \code{FALSE}.
#'
#' @export
egf_has_random <- function(object) {
  UseMethod("egf_has_random", object)
}

#' @rdname egf_has_random
#' @export
egf_has_random.egf <- function(object) {
  egf_has_random(object$tmb_args$data)
}

#' @export
egf_has_random.tmb_data <- function(object) {
  ncol(object$Z) > 0L
}

#' Construct the combined model frame
#'
#' Joins in a single \link[=data.frame]{data frame} all mixed effects
#' \link[=model.frame]{model frames} and further variables originally
#' indicated by \code{\link{egf}} argument \code{append}.
#'
#' @param object An \code{"\link{egf}"} object.
#'
#' @details
#' If a variable name occurs in multiple mixed effects model frames,
#' then only one instance is retained. Except in unusual cases
#' (possible only if model formulae have different formula environments),
#' all instances of a variable name are identical, and no information is lost.
#'
#' Since the data frames being combined each correspond rowwise
#' to \code{object$frame_windows}, so does the result.
#'
#' @return
#' A \link[=data.frame]{data frame} combining
#' (in the sense of \code{\link{cbind}})
#' all data frames in \link{list} \code{object$frame_parameters}
#' and data frame \code{object$frame_append}.
#'
#' @export
egf_combine_frames <- function(object) {
  stopifnot(inherits(object, "egf"))
  l <- c(unname(object$frame_parameters), list(object$frame_append))
  res <- do.call(cbind, l)
  res[!duplicated(names(res))]
}
