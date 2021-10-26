#' Compute likelihood profiles
#'
#' Computes univariate likelihood profiles of fixed effects
#' coefficients, random effect covariance parameters, and
#' linear combinations thereof.
#'
#' @param fitted
#'   An \code{"\link{egf}"} object.
#' @param level
#'   A number in the interval (0,1) indicating a confidence level.
#'   Profiles will be computed up to a deviance
#'   of \code{\link{qchisq}(level, df = 1)}.
#' @param which
#'   An integer vector indexing \code{fitted$best[!fitted$random]}.
#'   Only indexed parameters are profiled.
#' @param A
#'   A numeric matrix with \code{sum(!fitted$random)} columns,
#'   Each row specifies a linear combination of the elements
#'   of \code{fitted$best[!fitted$random]} to be profiled.
#'   Ignored if \code{which} is non-\code{NULL}.
#' @param top
#'   A subset of \code{\link{egf_get_names_top}(fitted, link = TRUE)}
#'   naming top level nonlinear model parameters whose population
#'   fitted values (see \code{\link[=fitted.egf]{fitted}}) should be
#'   profiled.
#'   Ignored if \code{which} or \code{A} is non-\code{NULL}.
#' @param parallel
#'   An \code{"\link{egf_parallel}"} object defining options for \R level
#'   parallelization.
#' @param trace
#'   A logical flag.
#'   If \code{TRUE}, then basic tracing messages indicating progress
#'   are printed.
#'   Depending on \code{fitted$control$trace}, these may be mixed with
#'   optimizer output.
#' @param grid_len
#'   A positive integer. Step sizes chosen adaptively by
#'   \code{\link[TMB]{tmbprofile}} will generate approximately
#'   this many points on each side of a profile's minimum point.
#' @param subset
#'   An expression to be evaluated in
#'   \code{\link[=model.frame.egf]{model.frame}(fitted, "combined")}
#'   It must evaluate to to a valid index vector for the rows of the
#'   data frame and, in turn, fitting windows.
#'   Only population fitted values for indexed windows are profiled.
#'   The default (\code{NULL}) is to consider all windows.
#'   Ignored if \code{which} or \code{A} is non-\code{NULL}.
#' @param append
#'   An expression indicating variables in
#'   \code{\link[=model.frame.egf]{model.frame}(object, "combined")}
#'   to be included with the result.
#'   The default (\code{NULL}) is to append nothing.
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
#' See topic \code{\link{egf_eval}} for details on nonstandard evaluation
#' of \code{subset} and \code{append}.
#'
#' @return
#' A data frame inheriting from class \code{"egf_profile"}, with variables:
#' \item{top}{
#'   (\code{top}-based calls only.)
#'   Top level nonlinear model parameter,
#'   from \code{\link{egf_get_names_top}(fitted, link = link)}.
#' }
#' \item{ts}{
#'   (\code{top}-based calls only.)
#'   Time series, from
#'   \code{levels(\link[=model.frame.egf]{model.frame}(fitted)$ts)}.
#' }
#' \item{window}{
#'   (\code{top}-based calls only.)
#'   Fitting window, from
#'   \code{levels(\link[=model.frame.egf]{model.frame}(fitted)$window)}.
#' }
#' \item{linear_combination}{
#'   Row index of linear combination, from \code{seq_len(nrow(A))}.
#' }
#' \item{value}{
#'   Value of linear combination being profiled.
#' }
#' \item{deviance}{
#'   Deviance of the restricted model that assumes \code{value}
#'   for the linear combination being profiled.
#' }
#' \code{A}, \code{x = fitted$best[!fitted$random]}, and \code{level}
#' are retained as \link{attributes}.
#'
#' @examples
#' fitted <- egf_cache("egf-1.rds")
#' zz <- egf_cache("profile-egf-1.rds", {
#'   profile(fitted, subset = (country == "A" & wave == 1))
#' })
#' str(zz)
#'
#' @seealso \code{\link{confint.egf_profile}}, \code{\link{plot.egf_profile}}
#' @export
#' @importFrom Matrix sparseMatrix
#' @importMethodsFrom Matrix t tcrossprod rowSums diag
#' @importFrom methods is
#' @importFrom stats vcov model.frame
#' @import parallel
profile.egf <- function(fitted,
                        level = 0.95,
                        which = NULL,
                        A = NULL,
                        top = egf_get_names_top(fitted, link = TRUE),
                        parallel = egf_parallel(),
                        trace = FALSE,
                        grid_len = 12,
                        subset = NULL,
                        append = NULL,
                        .subset = NULL,
                        .append = NULL,
                        ...) {
  stopifnot(
    is_number_in_interval(level, 0, 1, "()"),
    inherits(parallel, "egf_parallel"),
    is_true_or_false(trace),
    is_number_in_interval(grid_len, 1, Inf, "[)")
  )
  n <- sum(!fitted$random)

  ## If profiling user-specified elements of 'c(beta, theta)'
  if (!is.null(which)) {
    method <- "which"
    eval(bquote(stopifnot(
      is.numeric(which),
      which %in% seq_len(.(n))
    )))
    which <- unique(which)
    A <- sparseMatrix(i = seq_len(m), j = which, x = 1, dims = c(length(which), n))

  ## If profiling user-specified linear combinations
  ## of elements of 'c(beta, theta)'
  } else if (!is.null(A)) {
    method <- "A"
    if (!is(A, "dMatrix")) {
      stopifnot(is.numeric(A))
      if (is.null(dim(A))) {
        dim(A) <- c(1L, length(A))
      } else {
        stopifnot(is.matrix(A))
      }
    }
    eval(bquote(stopifnot(
      nrow(A) > 0L,
      ncol(A) == .(n),
      is.finite(A),
      rowSums(abs(A)) > 0
    )))

  ## If profiling population fitted values
  ## of top level nonlinear model parameters
  } else if (!is.null(top)) {
    method <- "top"

    names_top <- egf_get_names_top(fitted, link = TRUE)
    top <- unique(match.arg(top, names_top, several.ok = TRUE))

    frame_windows <- model.frame(fitted, "windows")
    frame_combined <- model.frame(fitted, "combined")
    subset <- if (is.null(.subset)) substitute(subset) else .subset
    subset <- egf_eval_subset(subset, frame_combined, parent.frame())
    append <- if (is.null(.append)) substitute(append) else .append
    append <- egf_eval_append(append, frame_combined, baseenv())

    l <- egf_preprofile(fitted, subset = subset, top = top)
    Y <- l$Y
    A <- l$A

  ## Otherwise
  } else {
    stop("One of 'which', 'A', and 'top' must be non-NULL.")
  }

  ## Covariance matrix of 'c(beta, theta)'
  V <- unclass(vcov(fitted))
  m <- nrow(A)
  if (method == "which") {
    a <- which
    h <- 0.25 * sqrt(diag(V)[which])
    s <- "name"
  } else {
    ## Covariance matrix of 'A %*% c(beta, theta)'
    V <- A %*% tcrossprod(V, A)
    a <- lapply(seq_len(m), function(i) A[i, ])
    h <- 0.25 * sqrt(diag(V))
    s <- "lincomb"
  }

  ytol <- 0.5 * qchisq(level, df = 1) # y := nll_restricted - nll_minimum = 0.5 * deviance
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
    args <- egf_tmb_remake_args(fitted$tmb_out, par = fitted$best)

    ## Retrieve path to shared object for loading
    dll <- system.file("libs", TMB::dynlib("epigrowthfit"), package = "epigrowthfit", mustWork = TRUE)

    cl <- parallel$cl
    if (is.null(cl)) {
      cl <- do.call(makePSOCKcluster, parallel$args)
      on.exit(stopCluster(cl), add = TRUE)
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
    if (given_outfile <- nzchar(parallel$outfile)) {
      outfile <- file(parallel$outfile, open = "wt")
      sink(outfile, type = "output")
      sink(outfile, type = "message")
      on.exit(add = TRUE, {
        sink(type = "message")
        sink(type = "output")
      })
    }
    if ((onomp <- TMB::openmp(n = NULL)) > 0L) {
      TMB::openmp(n = nomp)
      on.exit(TMB::openmp(n = onomp), add = TRUE)
    }
    obj <- fitted$tmb_out
    res <- switch(parallel$method,
      multicore = do.call(mcMap, c(list(f = do_profile, i = seq_len(m), a = a, h = h), parallel$args)),
      serial = Map(do_profile, i = seq_len(m), a = a, h = h)
    )
  }

  nr <- vapply(res, nrow, 0L)
  res <- data.frame(
    linear_combination = rep.int(gl(m, 1L), nr),
    do.call(rbind, res),
    row.names = NULL
  )
  if (method == "top") {
    i <- rep.int(rep.int(subset, length(top)), nr)
    res <- data.frame(
      top = rep.int(rep(factor(top, levels = names_top), each = length(subset)), nr),
      frame_windows[i, c("ts", "window"), drop = FALSE],
      res,
      frame_combined[i, append, drop = FALSE],
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    res$value <- res$value + rep.int(as.double(Y), nr)
  }
  attr(res, "A") <- A
  attr(res, "x") <- fitted$best[!fitted$random]
  attr(res, "level") <- level
  attr(res, "method") <- method
  class(res) <- c("egf_profile", oldClass(res))
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
#'   A logical flag. If \code{FALSE} and \code{object} supplies
#'   likelihood profiles of population fitted values of top level
#'   nonlinear model parameters, then confidence intervals
#'   on inverse link-transformed fitted values are returned.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' Each supplied likelihood profile
#' (level of factor \code{object$linear_combination}),
#' is linearly interpolated to approximate the two solutions
#' of \code{deviance(value) = \link{qchisq}(level, df = 1)}.
#' These provide the lower and upper confidence limits of interest
#' (see \href{https://en.wikipedia.org/wiki/Wilks'_theorem}{Wilks' theorem}).
#'
#' @return
#' A data frame with one row per supplied profile, and variables:
#' \item{linear_combination}{
#'   Row index of linear combination that was profiled,
#'   from \code{seq_len(nrow(attr(object, "A")))}.
#' }
#' \item{estimate, lower, upper}{
#'   Estimate of linear combination and approximate lower
#'   and upper confidence limits, inverse-link transformed
#'   if \code{link = FALSE}.
#' }
#' \code{level} is retained as an attribute of the result.
#' So are attributes \code{A} and \code{x} of \code{object}.
#'
#' @examples
#' object <- egf_cache("profile-egf-1.rds")
#' zz <- confint(object)
#' str(zz)
#'
#' @export
#' @importFrom stats qchisq approx
confint.egf_profile <- function(object, parm, level = attr(object, "level"), link = TRUE, ...) {
  stopifnot(
    is_true_or_false(link),
    is_number_in_interval(level, 0, attr(object, "level"), "(]")
  )
  q <- qchisq(level, df = 1)
  method <- attr(object, "method")

  do_solve <- function(d) {
    i_min <- which.min(d$deviance)
    i_left <- seq_len(i_min)
    i_right <- seq.int(i_min, nrow(d))
    estimate <- d$value[i_min]
    lower <- approx(x = d$deviance[i_left],  y = d$value[i_left],  xout = q)$y
    upper <- approx(x = d$deviance[i_right], y = d$value[i_right], xout = q)$y
    d1 <- d[1L, , drop = FALSE]
    m <- match(c("value", "deviance"), names(d1), 0L)
    data.frame(
      d1[seq_len(m[1L] - 1L)],
      estimate,
      lower,
      upper,
      d1[-seq_len(m[2L])],
      row.names = FALSE,
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
#' @param level
#'   A numeric vector with elements in (0,1). If \code{sqrt = FALSE},
#'   then line segments are drawn to show the intersection of the profile
#'   with lines at \code{deviance = \link{qchisq}(level, df = 1)}.
#' @param sqrt
#'   A logical flag.
#'   If \code{TRUE}, then square root-transformed deviance is plotted.
#' @param subset
#'   An expression to be evaluated in \code{x}.
#'   It must evaluate to a valid index vector for the rows of \code{x}.
#'   Only indexed profiles are plotted.
#'   The default (\code{NULL}) is to plot all profiles.
#' @param ...
#'   Graphical parameters passed to \code{\link{plot}}.
#'
#' @details
#' See topic \code{\link{egf_eval}} for details on nonstandard evaluation
#' of \code{subset}.
#'
#' @return
#' \code{NULL} (invisibly).
#'
#' @examples
#' x <- egf_cache("profile-egf-1.rds")
#' plot(x, type = "o", bty = "u", las = 1, main = "")
#'
#' @export
#' @import graphics
#' @importFrom stats confint
plot.egf_profile <- function(x, level = attr(x, "level"), sqrt = FALSE, subset = NULL, ...) {
  subset <- egf_eval_subset(substitute(subset), x, parent.frame())
  subset <- match(levels(factor(x$linear_combination[subset])), levels(x$linear_combination))

  stopifnot(is_true_or_false(sqrt))
  f <- if (sqrt) base::sqrt else identity
  do_segments <-
    !sqrt &&
    is.numeric(level) &&
    length(level) > 0L &&
    !all(is.na(level))
  if (do_segments) {
    level <- level[!is.na(level)]
    stopifnot(level > 0, level <= attr(x, "level"))
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
      main <- character(length(subset))
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
      segments(
        x0 = v_lower[[i]],
        x1 = v_upper[[i]],
        y0 = h,
        y1 = h,
        lty = 3
      )
      segments(
        x0 = c(v_lower[[i]], v_upper[[i]]),
        x1 = c(v_lower[[i]], v_upper[[i]]),
        y0 = usr[3L],
        y1 = rep.int(h, 2L),
        lty = 3
      )
      text(
        x = mean(usr[1:2]),
        y = h,
        labels = sprintf("%.3g%%", 100 * level),
        pos = 3,
        offset = 0.1
      )
    }
  }

  invisible(NULL)
}
