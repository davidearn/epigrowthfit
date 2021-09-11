#' Compute likelihood profiles
#'
#' Computes univariate likelihood profiles of fixed effects
#' coefficients, random effects covariance parameters, and
#' linear combinations thereof.
#'
#' @param fitted
#'   An \code{"\link{egf}"} object.
#' @param which
#'   An integer vector indexing coefficients in \code{fitted$best}
#'   to be profiled. Must be a subset of \code{fitted$nonrandom}.
#' @param A
#'   A \link{numeric} \link{matrix} with
#'   \code{\link{length}(fitted$nonrandom)} columns.
#'   Each row specifies a linear combination of the elements
#'   of \code{fitted$best[fitted$nonrandom]} to be profiled.
#'   Ignored if \code{which} is non-\code{\link{NULL}}.
#' @param top
#'   A subset of \code{\link{egf_get_names_top}(fitted, link = TRUE)}
#'   naming top level nonlinear model parameters whose population
#'   fitted values (see \code{\link{fitted.egf}}) should be profiled.
#'   Ignored if \code{which} or \code{A} is non-\code{\link{NULL}}.
#' @param subset
#'   An expression to be evaluated in the combined model frame
#'   (see \code{\link{egf_combine_frames}}). It must evaluate to
#'   to a valid index vector for the rows of the data frame
#'   (see \code{\link{[.data.frame}}), and thus fitting windows.
#'   Only population fitted values for indexed windows are profiled.
#'   The default (\code{\link{NULL}}) is to consider all windows.
#'   Ignored if \code{which} or \code{A} is non-\code{\link{NULL}}.
#' @param append
#'   An expression indicating variables in the combined model frame
#'   (see \code{\link{egf_combine_frames}}) to be included with the
#'   result. The default (\code{\link{NULL}}) is to append nothing.
#'   Ignored if \code{which} or \code{A} is non-\code{\link{NULL}}.
#' @param level_max
#'   A number in the interval (0,1) indicating a confidence level.
#'   Profiles will be computed up to a deviance
#'   of \code{\link{qchisq}(level_max, df = 1)}.
#' @param grid_len
#'   A positive integer. Step sizes chosen adaptively by
#'   \code{\link[TMB]{tmbprofile}} will generate approximately
#'   this many points on each side of a profile's minimum point.
#' @param trace
#'   A \link{logical} flag.
#'   If \code{TRUE}, then basic tracing messages indicating progress
#'   are printed.
#'   Depending on \code{fitted$control$trace}, these may be mixed with
#'   optimizer output.
#' @param parallel
#'   An \code{"\link{egf_parallel}"} object defining options for \R level
#'   parallelization.
#' @param .subset,.append
#'   Index vectors to be used (if non-\code{\link{NULL}}) in place of
#'   the result of evaluating \code{subset} and \code{append}.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' Computation of likelihood profiles is expensive as it requires estimation
#' of many restricted models. It is parallelized at the C++ level when there
#' is OpenMP support and \code{fitted$control$omp_num_threads} is set to
#' an integer greater than 1. If there is no OpenMP support, then computation
#' can still be parallelized at the \R level with appropriate setting of
#' \code{parallel}.
#'
#' \code{which} is mapped to an \code{A} matrix composed of unit row vectors,
#' with 1 at array index \code{[i, which[i]]} for all \code{i}, such that
#' \code{A \link[=matmult]{\%*\%} fitted$best[fitted$nonrandom]} is precisely
#' \code{fitted$best[which]}.
#'
#' \code{top} and \code{subset} are mapped to a block \code{A} matrix.
#' Block \code{k} (from the top) is composed of rows (those indexed
#' by \code{subset}) of the fixed effects design matrix for parameter
#' \code{top[k]} .
#' \code{A \link[=matmult]{\%*\%} fitted$best[fitted$nonrandom]}
#' is a vector listing the population fitted values for each parameter
#' named in \code{top}, for each fitting window indexed by \code{subset}.
#'
#' See topic \code{\link{egf_eval}} for details on nonstandard evaluation
#' of \code{subset} and \code{append}.
#'
#' @return
#' A \link[=data.frame]{data frame} inheriting from \link{class}
#' \code{"egf_profile"}, with variables:
#' \item{top}{
#'   (\code{top}-based calls only.)
#'   Top level nonlinear model parameter,
#'   from \code{\link{egf_get_names_top}(fitted, link = link)}.
#' }
#' \item{ts}{
#'   (\code{top}-based calls only.)
#'   Time series, from \code{\link{levels}(fitted$frame_windows$ts)}.
#' }
#' \item{window}{
#'   (\code{top}-based calls only.)
#'   Fitting window, from \code{\link{levels}(fitted$frame_windows$window)}.
#' }
#' \item{linear_combination}{
#'   Row index of linear combination,
#'   from \code{\link{seq_len}(\link{nrow}(A))}.
#' }
#' \item{value}{
#'   Value of linear combination being profiled.
#' }
#' \item{deviance}{
#'   Deviance of the restricted model that assumes \code{value}
#'   for the linear combination being profiled.
#' }
#' \code{A}, \code{x = fitted$best[fitted$nonrandom]}, and \code{level_max}
#' are retained as \link{attributes}.
#'
#' @examples
#' example("egf", "epigrowthfit")
#' zz <- profile(object, subset = (country == "A" & wave == 1))
#' str(zz)
#'
#' @seealso \code{\link{confint.egf_profile}}, \code{\link{plot.egf_profile}}
#' @export
#' @importFrom Matrix sparseMatrix KhatriRao
#' @importMethodsFrom Matrix diag
#' @importFrom methods as is
#' @importFrom stats vcov
#' @import parallel
profile.egf <- function(fitted,
                        which = NULL,
                        A = NULL,
                        top = egf_get_names_top(fitted, link = TRUE),
                        subset = NULL,
                        append = NULL,
                        level_max = 0.95,
                        grid_len = 12,
                        trace = TRUE,
                        parallel = egf_parallel(),
                        .subset = NULL,
                        .append = NULL,
                        ...) {
  stop_if_not_number_in_interval(level_max, 0, 1, "()")
  stop_if_not_number_in_interval(grid_len, 1, Inf, "[)")
  stop_if_not_true_false(trace)
  stopifnot(inherits(parallel, "egf_parallel"))
  n <- length(fitted$nonrandom)

  ## If profiling user-specified elements of 'c(beta, theta)'
  if (!is.null(which)) {
    method <- "which"
    stopifnot(
      is.numeric(which),
      which %in% fitted$nonrandom
    )
    which <- unique(which)
    m <- length(which)
    A <- sparseMatrix(i = seq_len(m), j = which, x = 1, dims = c(m, n))

  ## If profiling user-specified linear combinations
  ## of elements of 'c(beta, theta)'
  } else if (!is.null(A)) {
    method <- "A"
    stopifnot(is.numeric(A))
    if (is.null(dim(A))) {
      dim(A) <- c(1L, length(A))
    } else {
      stopifnot(is.matrix(A) || is(A, "dMatrix"))
    }
    stopifnot(
      ncol(A) == length(fitted$nonrandom),
      is.finite(A),
      rowSums(abs(A)) > 0
    )
    m <- nrow(A)

  ## If profiling population fitted values
  ## of nonlinear and dispersions model parameters
  } else if (!is.null(top)) {
    method <- "top"
    names_top <- egf_get_names_top(fitted, link = TRUE)
    top <- unique(match.arg(top, names_top, several.ok = TRUE))
    combined <- egf_combine_frames(fitted)
    subset <- if (is.null(.subset)) substitute(subset) else .subset
    subset <- egf_eval_subset(subset, combined, parent.frame())
    append <- if (is.null(.append)) substitute(append) else .append
    append <- egf_eval_append(append, combined, baseenv())

    p <- length(top)
    N <- length(subset)
    m <- p * N

    f <- factor(fitted$info$X$top, levels = top)
    J <- as(f, "sparseMatrix")
    X <- fitted$tmb_out$env$data$X[subset, , drop = FALSE]
    A <- KhatriRao(J, X)
    if (egf_has_random(fitted)) {
      index <- grepl("^beta\\[", names(fitted$best)[fitted$nonrandom])
      A@p <- c(0L, replace(rep_len(NA_integer_, n), index, A@p[-1L]))
      A@p <- locf(A@p)
      A@Dim <- c(m, n)
    }

  ## Otherwise
  } else {
    stop("One of 'A', 'which', and 'top' must be non-NULL.")
  }

  ## Covariance matrix of 'c(beta, theta)'
  V <- vcov(fitted, full = TRUE)
  if (method == "which") {
    a <- which
    h <- sqrt(diag(V)[which]) / 4
    s <- "name"
  } else {
    ## Covariance matrix of 'A %*% c(beta, theta)'
    V <- A %*% unclass(V) %*% t(A)
    a <- lapply(seq_len(nrow(A)), function(i) A[i, ])
    h <- sqrt(diag(V)) / 4
    s <- "lincomb"
  }

  ytol <- 0.5 * qchisq(level_max, df = 1) # y := nll_restricted - nll_minimum = 0.5 * deviance
  ystep <- ytol / grid_len
  nomp <- fitted$control$omp_num_threads

  do_profile <- function(i, a, h) {
    if (trace) {
      cat(sprintf("Computing likelihood profile %d of %d...\n", i, m))
    }
    args <- list(obj = obj, h = h, ytol = ytol, ystep = ystep)
    args[[s]] <- a
    res <- do.call(TMB::tmbprofile, args)
    i_min <- which.min(res[[2L]])
    res[[2L]] <- 2 * (res[[2L]] - res[i_min, 2L]) # deviance = 2 * (nll_restricted - nll_minimum)
    names(res) <- c("value", "deviance")
    res[-i_min, , drop = FALSE] # 'tmbprofile' duplicates this row
  }

  if (parallel$method == "snow") {
    environment(do_profile) <- .GlobalEnv

    ## Reconstruct list of arguments to 'MakeADFun' from object internals
    ## for retaping
    args <- egf_tmb_remake_args(fitted)

    ## Retrieve path to shared object for loading
    dll <- system.file("libs", TMB::dynlib("epigrowthfit"), package = "epigrowthfit")

    if (is.null(parallel$cl)) {
      cl <- do.call(makePSOCKcluster, parallel$args)
      on.exit(stopCluster(cl), add = TRUE)
    } else {
      cl <- parallel$cl
    }
    clusterExport(cl,
      varlist = c("dll", "nomp", "args", "trace", "m", "s", "ytol", "ystep"),
      envir = environment()
    )
    clusterEvalQ(cl, {
      dyn.load(dll)
      if (TMB::openmp(n = NULL) > 0L) {
        TMB::openmp(n = nomp)
      }
      obj <- do.call(TMB::MakeADFun, args)
    })
    res <- clusterMap(cl, do_profile, i = seq_len(m), a = a, h = h)
  } else {
    obj <- fitted$tmb_out
    if ((onomp <- TMB::openmp(n = NULL)) > 0L) {
      TMB::openmp(n = nomp)
      on.exit(TMB::openmp(n = onomp), add = TRUE)
    }
    if (given_outfile <- nzchar(parallel$outfile)) {
      outfile <- file(parallel$outfile, open = "wt")
      sink(outfile, type = "output")
      sink(outfile, type = "message")
    }
    res <- tryCatch(
      switch(parallel$method,
        multicore = do.call(mcMap, c(list(f = do_profile, i = seq_len(m), a = a, h = h), parallel$args)),
        serial = Map(do_profile, i = seq_len(m), a = a, h = h)
      ),
      error = function(cond) {
        if (given_outfile) {
          sink(type = "message")
          sink(type = "output")
        }
        stop(cond)
      }
    )
  }

  nr <- vapply(res, nrow, 0L)
  res <- data.frame(
    linear_combination = rep.int(gl(m, 1L), nr),
    do.call(rbind, res),
    row.names = NULL
  )
  if (method == "top") {
    i <- rep.int(rep.int(subset, p), nr)
    res <- data.frame(
      top = rep.int(rep.int(factor(top, levels = names_top), N), nr),
      fitted$frame_windows[i, c("ts", "window"), drop = FALSE],
      res,
      combined[i, append, drop = FALSE],
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  attr(res, "A") <- A
  attr(res, "x") <- fitted$best[fitted$nonrandom]
  attr(res, "level_max") <- level_max
  attr(res, "method") <- method
  class(res) <- c("egf_profile", "data.frame")
  res
}

#' Confidence intervals from likelihood profiles
#'
#' Computes confidence intervals on fixed effects coefficients,
#' random effects covariance parameters, and linear combinations
#' thereof from their univariate likelihood profiles.
#'
#' @param object
#'   An \code{"\link[=profile.egf]{egf_profile}"} object.
#' @param parm
#'   Unused argument included for generic consistency.
#' @param level
#'   A number in the interval (0,1) indicating a confidence level.
#' @param link
#'   A \link{logical} flag. If \code{FALSE} and \code{object}
#'   supplies likelihood profiles of population fitted values
#'   of top level nonlinear model parameters, then confidence
#'   intervals on inverse link-transformed fitted values are
#'   returned.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' Each supplied likelihood profile
#' (level of \link{factor} \code{object$linear_combination}),
#' is linearly interpolated to approximate the two solutions
#' of \code{deviance(value) = \link{qchisq}(level, df = 1)}.
#' These provide the lower and upper confidence limits of interest
#' (see \href{https://en.wikipedia.org/wiki/Wilks'_theorem}{Wilks' theorem}).
#'
#' @return
#' A \link[=data.frame]{data frame} with one row per supplied profile,
#' and variables:
#' \item{linear_combination}{
#'   Row index of linear combination that was profiled,
#'   from \code{\link{seq_len}(\link{nrow}(\link{attr}(object, "A")))}.
#' }
#' \item{estimate, lower, upper}{
#'   Estimate of linear combination and approximate lower
#'   and upper confidence limits, inverse-link transformed
#'   if \code{link = FALSE}.
#' }
#' \code{level} is retained as an \link[=attributes]{attribute} of the result.
#' So are attributes \code{A} and \code{x} of \code{object}.
#'
#' @examples
#' example("profile.egf", "epigrowthfit")
#' confint(zz, link = TRUE)
#' confint(zz, link = FALSE)
#'
#' @export
#' @importFrom stats qchisq approx
confint.egf_profile <- function(object, parm, level = attr(object, "level_max"), link = TRUE, ...) {
  stop_if_not_true_false(link)
  stop_if_not_number_in_interval(level, 0, 1, "()")
  stopifnot(level <= attr(object, "level_max"))
  q <- qchisq(level, df = 1)
  method <- attr(object, "method")

  s <- c("top", "ts", "window", "linear_combination", "value", "deviance")
  if (method == "top") {
    j1 <- 1:4
    j2 <- -match(s, names(object), 0L)
  } else {
    j1 <- 4L
    j2 <- -match(s[4:6], names(object), 0L)
  }
  do_solve <- function(d) {
    i_min <- which.min(d$deviance)
    i_left <- seq_len(i_min)
    i_right <- seq.int(i_min, nrow(d))
    estimate <- d$value[i_min]
    lower <- approx(x = d$deviance[i_left],  y = d$value[i_left],  xout = q)$y
    upper <- approx(x = d$deviance[i_right], y = d$value[i_right], xout = q)$y
    data.frame(
      d[1L, j1, drop = FALSE],
      estimate,
      lower,
      upper,
      d[1L, j2, drop = FALSE],
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  res <- do.call(rbind, by(object, object$linear_combination, do_solve, simplify = FALSE))
  if (method == "top" && !link) {
    elu <- c("estimate", "lower", "upper")
    res[elu] <- in_place_ragged_apply(res[elu], res$top,
      f = lapply(egf_link_extract(levels(res$top)), egf_link_match, inverse = TRUE)
    )
    levels(res$top) <- egf_link_remove(levels(res$top))
  }
  res$linear_combination <- as.integer(as.character(res$linear_combination))
  row.names(res) <- NULL
  attr(res, "A") <- attr(object, "A")
  attr(res, "x") <- attr(object, "x")
  attr(res, "level") <- level
  res
}

#' Plot likelihood profiles
#'
#' A method for plotting likelihood profiles.
#'
#' @param x
#'   An \code{"\link[=profile.egf]{egf_profile}"} object.
#' @param subset
#'   An expression to be evaluated in \code{x}.
#'   It must evaluate to a valid index vector for the rows of \code{x}
#'   (see \code{\link{[.data.frame}}).
#'   Only indexed profiles are plotted.
#'   The default (\code{\link{NULL}}) is to plot all profiles.
#' @param sqrt
#'   A \link{logical} flag.
#'   If \code{TRUE}, then square root-transformed deviance is plotted.
#' @param level
#'   A \link{numeric} vector with elements in (0,1). If \code{sqrt = FALSE},
#'   then line segments are drawn to show the intersection of the profile
#'   with lines at \code{deviance = \link{qchisq}(level, df = 1)}.
#' @param .segments,.text
#'   \link[=list]{List}s of optional graphical parameters passed to
#'   \code{\link{segments}} and \code{\link{text}} when line segments
#'   are drawn.
#' @param ...
#'   Optional graphical parameters passed to \code{\link{plot}}.
#'
#' @details
#' See topic \code{\link{egf_eval}} for details on nonstandard evaluation
#' of \code{subset}.
#'
#' @return
#' \code{\link{NULL}} (invisibly).
#'
#' @examples
#' example("profile.egf", "epigrowthfit")
#' plot(zz, type = "o", bty = "u", las = 1, main = "", .segments = list(lty = 3))
#'
#' @export
#' @import graphics
#' @importFrom stats confint
plot.egf_profile <- function(x, subset = NULL, sqrt = FALSE,
                             level = attr(x, "level_max"),
                             .segments = list(), .text = list(), ...) {
  subset <- egf_eval_subset(substitute(subset), x, parent.frame())
  subset <- match(levels(factor(x$linear_combination[subset])), levels(x$linear_combination))

  stop_if_not_true_false(sqrt)
  f <- if (sqrt) base::sqrt else identity
  do_segments <-
    !sqrt &&
    is.numeric(level) &&
    length(level) > 0L &&
    !all(is.na(level))
  if (do_segments) {
    level <- level[!is.na(level)]
    stopifnot(
      level > 0,
      level <= attr(x, "level_max"),
      is.list(.segments),
      is.list(.text)
    )
    ## Line segments at heights 'h' in all plots
    h <- qchisq(level, df = 1)
    ## Line segment 'j' to start at 'v_lower[[i]][j]'
    ## and end at 'v_upper[[i]][j]' in plot 'i'
    ci <- lapply(level, function(p) confint(x, level = p))
    v_lower <- lapply(subset, function(i) vapply(ci, `[`, 0, i, "lower"))
    v_upper <- lapply(subset, function(i) vapply(ci, `[`, 0, i, "upper"))
  }

  dots <- list(...)
  method <- attr(x, "method")
  if (is.null(dots[["main"]])) {
    if (method == "top") {
      main <- c(tapply(as.character(x$window), x$linear_combination, `[[`, 1L))[subset]
    } else {
      main <- rep_len("", length(subset))
    }
  } else {
    main <- rep_len(dots$main, length(subset))
  }
  if (is.null(dots[["xlab"]])) {
    if (method == "top") {
      xlab <- c(tapply(as.character(x$top), x$linear_combination, `[[`, 1L))[subset]
    } else {
      xlab <- paste("linear combination", levels(x$linear_combination)[subset])
    }
  } else {
    xlab <- rep_len(dots$xlab, length(subset))
  }
  dots$main <- dots$xlab <- NULL
  if (is.null(dots[["ylab"]])) {
    dots$ylab <- "deviance"
    if (sqrt) {
      dots$ylab <- as.expression(bquote(sqrt(.(dots$ylab))))
    }
  }
  if (is.null(dots[["ylim"]])) {
    dots$ylim <- c(0, f(max(x$deviance, na.rm = TRUE)))
  }

  for (i in seq_along(subset)) {
    args <- list(
      formula = f(deviance) ~ value,
      data = x,
      subset = (unclass(x$linear_combination) == subset[i]),
      main = main[i],
      xlab = xlab[i]
    )
    do.call(plot, c(args, dots))
    if (do_segments) {
      usr <- par("usr")
      args <- list(
        x0 = v_lower[[i]],
        x1 = v_upper[[i]],
        y0 = h,
        y1 = h
      )
      do.call(segments, c(args, .segments))
      args <- list(
        x0 = c(v_lower[[i]], v_upper[[i]]),
        x1 = c(v_lower[[i]], v_upper[[i]]),
        y0 = usr[3L],
        y1 = rep.int(h, 2L)
      )
      do.call(segments, c(args, .segments))
      args <- list(
        x = mean(usr[1:2]),
        y = h,
        labels = sprintf("%.3g%%", 100 * level),
        pos = 3,
        offset = 0.1
      )
      do.call(text, c(args, .text))
    }
  }
  invisible(NULL)
}
