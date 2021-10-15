#' Define an epidemic model
#'
#' Sets flags defining the top level nonlinear model of epidemic growth
#' to be estimated by \code{\link{egf}}.
#'
#' @param curve
#'   A character string specifying a model for expected cumulative
#'   disease incidence as a function of time.
#' @param excess
#'   A logical flag. If \code{TRUE}, then a constant baseline mortality rate
#'   is estimated. Set to \code{TRUE} if what is observed is multiple causes
#'   mortality rather than disease mortality or disease incidence.
#' @param family
#'   A character string specifying a family of discrete probability
#'   distributions assigned to observations of disease incidence.
#' @param day_of_week
#'   An integer flag. If positive, then day of week effects are estimated
#'   as offsets relative to the indicated day of week
#'   (Sunday if \code{day_of_week = 1}, Monday if \code{day_of_week = 2},
#'   and so on). Logical values are coerced to integer.
#'
#' @return
#' A list inheriting from class \code{"egf_model"} containing the arguments
#' (after possible matching and coercion).
#'
#' @examples
#' model <- egf_model()
#' str(model)
#'
#' @family behaviour-defining functions
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
  class(res) <- "egf_model"
  res
}

#' Define control parameters
#'
#' Set parameters controlling the behaviour of \code{\link{egf}}.
#'
#' @param optimizer
#'   An \code{"\link{egf_optimizer}"} object, specifying an "outer"
#'   optimization method.
#' @param inner_optimizer
#'   An \code{"\link{egf_inner_optimizer}"} object, specifying an "inner"
#'   optimization method, or a list of such objects, in which case the
#'   listed methods are tried in order until one succeeds. (If none succeed,
#'   then a warning is issued.)
#' @param trace
#'   An integer flag determining the amount of tracing performed
#'   (see Details). Logical values are coerced to integer.
#' @param profile
#'   A logical flag. If \code{TRUE}, then fixed effect parameters are profiled
#'   out of the likelihood, which may stabilize optimization for models
#'   with many fixed effects. This feature should be considered experimental,
#'   and in fact may \emph{de}stabilize optimization, as it may rely on
#'   assumptions about the optimization problem that are not necessarily
#'   satisfied by the nonlinear mixed effects models fit by \code{\link{egf}}.
#' @param sparse_X
#'   A logical flag. If \code{TRUE}, then the fixed effects
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
#' @section Warning:
#' Setting \code{trace > 0L} and \code{omp_num_threads > 0L}
#' simultaneously should be considered dangerous on builds of
#' \pkg{epigrowthfit} obtained from CRAN. These builds print
#' using R API that is not thread-safe. R API is avoided on
#' builds from source.
#'
#' @return
#' A list inheriting from class \code{"egf_control"} containing the arguments
#' (after possible coercion).
#'
#' @examples
#' control <- egf_control()
#' str(control)
#'
#' @family behaviour-defining functions
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
  class(res) <- "egf_control"
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
#'   A function performing optimization. The outer optimization permits
#'   \code{\link{optim}}, \code{\link{nlminb}}, and \code{\link{nlm}}
#'   and any \code{optim}-like function. An \code{optim}-like function
#'   is a function \code{f} such that:
#'   (i) the first three arguments of \code{f} specify an initial parameter
#'   vector, an objective function, and a gradient function, respectively;
#'   (ii) \code{f} accepts \code{control} as a fourth (or later) argument;
#'   and
#'   (iii) \code{f} returns a list with elements \code{par}, \code{value},
#'   \code{convergence}, and \code{message}.
#'   The inner optimization permits \code{optim} and \code{\link[TMB]{newton}}
#'   only.
#' @param args
#'   A list of arguments to \code{f} other than \code{control}.
#'   If \code{f = \link{optim}} and \code{args} does not have \code{method}
#'   as an element, then \code{method = "BFGS"} is appended.
#' @param control
#'   A list of control parameters to be assigned to argument \code{control}
#'   of \code{f}.
#'
#' @return
#' \code{egf_optimizer} returns a list inheriting from class
#' \code{"egf_optimizer"}, with elements:
#' \item{f}{
#'   An \code{\link{optim}}-like \link{function}, typically the result
#'   of wrapping the supplied optimizer to make it \code{optim}-like.
#' }
#' \item{args}{
#'   The supplied list of arguments
#'   (after possible deletion of elements with reserved names).
#' }
#' \item{control}{
#'   The supplied list of control parameters.
#' }
#'
#' \code{egf_inner_optimizer} returns a list inheriting from class
#' \code{"egf_inner_optimizer"}, with elements:
#' \item{method}{
#'   A character string. This is \code{args$method} if \code{f = \link{optim}}
#'   and \code{"newton"} if \code{f = \link[TMB]{newton}}.
#' }
#' \item{control}{
#'   A list. This is \code{control} if \code{f = \link{optim}} and \code{args}
#'   (after possible deletion of elements with reserved names)
#'   if \code{f = \link[TMB]{newton}}.
#'   To align the default behaviour of \code{newton} with that of \code{optim},
#'   \code{trace = 0} is set if not specified in \code{args}.
#' }
#'
#' @examples
#' optimizer <- egf_optimizer()
#' str(optimizer)
#'
#' inner_optimizer <- egf_inner_optimizer()
#' str(inner_optimizer)
#'
#' @name egf_optimizer
#' @seealso
#' \code{\link[TMB]{MakeADFun}} for some details about outer and inner
#' optimizations
NULL

#' @rdname egf_optimizer
#' @family behaviour-defining functions
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
    f_out <- tryCatch(
      eval(e),
      error = function(cond) {
        stop(wrap(
          "Unable to validate 'f' due to following error in test ", sQuote(deparse1(e)), ":\n\n",
          conditionMessage(cond)
        ))
      }
    )
    required <- c("par", "value", "convergence", "message")
    if (!(is.list(f_out) && all(required %in% names(f_out)))) {
      stop(wrap(
        "'f' must return a list with elements ",
        paste(sQuote(required), collapse = ", "),
        " but _does not_ for test ", sQuote(deparse1(e)), "."
      ))
    }
    f <- function(par, fn, gr, control, ...) {
      optimizer(par, fn, gr, control = control, ...)
    }
  }
  if (!is.null(names(args))) {
    reserved <- c("par", "fn", "gr", "control", "...", names(formals(optimizer))[1:3])
    args[reserved] <- NULL
  }

  res <- list(f = f, args = args, control = control)
  class(res) <- "egf_optimizer"
  res
}

#' @rdname egf_optimizer
#' @family behaviour-defining functions
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
  class(res) <- "egf_inner_optimizer"
  res
}

#' Define a parallelization method
#'
#' Defines instructions for parallelization by linking a method with options.
#'
#' @param method
#'   A character string indicating a method of parallelization.
#'   \code{"\link[=lapply]{serial}"} indicates no parallelization.
#'   \code{"\link[parallel:mclapply]{multicore}"} indicates \R level forking.
#'   It is intended for use from a terminal rather than a GUI.
#'   On Windows, \code{"multicore"} is equivalent to \code{"serial"}.
#'   \code{"\link[parallel:parLapply]{snow}"} indicates socket clusters.
#'   \code{"snow"} is supported on both Unix-alikes and Windows.
#' @param outfile
#'   A character string indicating a file path where console output
#'   should be diverted. An empty string indicates no diversion.
#'   If \code{method = "snow"}, then diversion may be necessary to view output.
#' @param cores
#'   A positive integer indicating a number of threads/processes
#'   to fork/spawn when \code{parallel != "serial"}. The maximum
#'   is typically \code{\link[parallel]{detectCores}(TRUE, FALSE)}.
#' @param args
#'   A list of optional arguments to
#'   \code{\link[parallel]{mclapply}} (\code{method = "multicore"}) or
#'   \code{\link[parallel]{makePSOCKcluster}} (\code{method = "snow"}).
#' @param cl
#'   An existing \link[parallel:makePSOCKcluster]{socket cluster}
#'   (\code{method = "snow"}).
#'   The default (\code{NULL}) is to create a new clusters as necessary
#'   and terminate them upon job completion.
#'   (If non-\code{NULL}, then \code{outfile}, \code{cores}, and \code{args}
#'   are ignored.)
#'
#' @details
#' For general information about parallelism in \R,
#' see \code{\link{vignette}("parallel", "parallel")}.
#'
#' @return
#' A list inheriting from class \code{"egf_parallel"}
#' containing the arguments (after possible matching and coercion).
#'
#' @examples
#' parallel <- egf_parallel()
#' str(parallel)
#'
#' @family behaviour-defining functions
#' @export
#' @importFrom TMB openmp
egf_parallel <- function(method = c("serial", "multicore", "snow"),
                         outfile = "",
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
  class(res) <- "egf_parallel"
  res
}

#' Define plot options
#'
#' Sets parameters controlling the graphical output of \code{\link{plot.egf}}.
#' Supplied values override package defaults
#' (retrievable as \code{defaults <- egf_plot_control()}),
#' which in turn override global defaults set via \code{\link{par}}.\cr\cr
#' Here, \code{x}, \code{type}, \code{time_as}, and \code{dt}
#' refer to the so-named arguments of \code{plot.egf}.
#'
#' @param window
#'   A named list of arguments to \code{\link{rect}}
#'   affecting the appearance of fitting windows.
#' @param data
#'   A named list of the form \code{list(main, short, long)}.
#'   \code{main} is a named list of arguments to \code{\link{points}}
#'   affecting the appearance of observed data.
#'   \code{short} and \code{long} are alternatives to \code{main}
#'   used for counts over intervals shorter or longer than \code{dt}
#'   when \code{type = "interval"}.
#'   \code{short} and \code{long} default to \code{main} (elementwise).
#' @param predict
#'   A named list of the form \code{list(estimate, ci)}.
#'   \code{estimate} and \code{ci} are named lists of arguments
#'   to \code{\link{lines}} and \code{\link{polygon}},
#'   affecting the appearance of predicted curves and corresponding
#'   confidence bands.
#'   \code{ci$col} defaults to \code{estimate$col} with added transparency.
#' @param asymptote
#'   A named list of arguments to \code{\link{segments}}
#'   affecting appearance of line segments drawn at
#'   \code{y = <initial exponential growth rate>}
#'   when \code{type = "rt"} and
#'   \code{x$model$curve = "logistic"} or \code{"richards"}.
#' @param box
#'   A named list of arguments to \code{\link{box}}
#'   affecting the appearance of the box drawn around the plot region
#'   on the device.
#' @param axis
#'   A named list of the form \code{\link{list}(x, y)}.
#'   \code{x} and \code{y} are named lists of arguments to \code{\link{axis}}
#'   affecting the appearance of the bottom and left axes.
#'   When \code{time_as = "Date"}, there are minor and major bottom axes.
#'   In this case, the major axis uses a modified version of \code{x}
#'   that tries to ensure that it is displayed below the minor axis in
#'   a slightly larger font.
#' @param title
#'   A named list of the form \code{list(main, sub, xlab, ylab, plab)}.
#'   The elements are named lists of arguments to \code{\link{title}}
#'   affecting the appearance of plot (sub)titles and axis and panel labels.
#'   \code{sub$adj} defaults to \code{main$adj}.
#' @param tdoubling
#'   A named list of the form \code{list(legend, estimate, ci)}.
#'   The elements are named lists of arguments to \code{\link{mtext}}
#'   affecting the appearance of initial doubling times printed in the
#'   top margin.
#' @param heat
#'   A named list of the form \code{list(pal, bg, ul, ips)}.
#'   \code{pal} is a named list of arguments to \code{\link{colorRamp}}
#'   defining a colour palette for heat maps.
#'   \code{bg} and \code{ul} are named lists of arguments to \code{\link{rect}}
#'   affecting the appearance of panel backgrounds and panel label underlays.
#'
#' @details
#' Setting an argument (or an element thereof, in the case of nested lists)
#' to \code{\link{NULL}} has the effect of suppressing the corresponding
#' plot element.
#'
#' For \code{type = "rt_heat"}, the only useful arguments are \code{axis},
#' \code{title}, and \code{heat}.
#'
#' @return
#' A named list.
#'
#' @export
egf_plot_control <- function(window, data, predict, asymptote,
                             box, axis, title, tdoubling, heat) {
  res <- list(
    window = list(
      col = alpha("#DDCC77", 0.25),
      border = NA
    ),
    data = list(
      main = list(
        pch = 21,
        col = "#BBBBBB",
        bg = "#DDDDDD"
      ),
      short = list(),
      long = list()
    ),
    predict = list(
      estimate = list(
        col = "#44AA99",
        lwd = 2
      ),
      ci = list(
        border = NA
      )
    ),
    asymptote = list(
      lty = "dotted",
      lwd = 2
    ),
    box = list(
      bty = "l"
    ),
    axis = list(
      x = list(
        gap.axis = 0
      ),
      y = list()
    ),
    title = list(
      main = list(adj = 0, xpd = NA),
      sub = list(mgp = c(0, 0, 0), xpd = NA),
      xlab = list(xpd = NA),
      ylab = list(xpd = NA),
      plab = list(col.lab = "white")
    ),
    tdoubling = list(
      legend = list(),
      estimate = list(),
      ci = list()
    ),
    heat = list(
      pal = list(
        colors = c(
          "#364B9A", "#4A7BB7", "#6EA6CD", "#98CAE1",
          "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366",
          "#F67E4B", "#DD3D2D", "#A50026"
        ),
        bias = 1,
        space = "rgb",
        interpolate = "linear"
      ),
      bg = list(
        col = "black",
        border = NA
      ),
      ul = list(
        col = alpha("black", 0.5),
        border = NA
      )
    )
  )

  ## Multi-assign from 'value' into 'default',
  ## recursively up to one level of nesting
  rmerge <- function(from, into, recursive = FALSE) {
    if (is.null(from)) {
      if (recursive) {
        into[] <- list(NULL)
        return(into)
      }
      return(NULL)
    }
    if (!is.list(from) || length(from) == 0L || is.null(names(from))) {
      return(into)
    }
    nel <- unique(names(from))
    if (recursive) {
      into[nel] <- Map(rmerge, from = from[nel], into = into[nel], recursive = FALSE)
    } else {
      into[nel] <- from[nel]
    }
    into
  }

  recursive <- c("data", "predict", "axis", "title", "tdoubling", "heat")
  nel <- names(match.call()[-1L])
  if (length(nel) > 0L) {
    res[nel] <- Map(rmerge,
      from = mget(nel, mode = "list", ifnotfound = list(list()), inherits = FALSE),
      into = res[nel],
      recursive = nel %in% recursive
    )
  }

  ## Some default values are conditional on supplied values
  for (s in c("short", "long")) {
    if (is.list(res$data[[s]])) {
      nel <- setdiff(names(res$data$main), names(res$data[[s]]))
      res$data[[s]][nel] <- res$data$main[nel]
    }
  }
  if (is.list(res$predict$ci) && is.null(res$predict$ci$col)) {
    res$predict$ci$col <- alpha(res$predict$estimate$col, 0.4)
  }
  adj <- res$title$main$adj
  if (is.numeric(adj) && length(adj) == 1L && is.finite(adj)) {
    if (is.list(res$title$sub) && is.null(res$title$sub$adj)) {
      res$title$sub$adj <- adj
    }
    if (is.list(res$tdoubling$legend) && is.null(res$tdoubling$legend$adj)) {
      res$tdoubling$legend$adj <- if (adj > 0.5) 0 else 1
    }
  }

  class(res) <- "egf_plot_control"
  res
}
