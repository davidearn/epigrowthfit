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
#' @param max_level
#'   A number in the interval (0,1) indicating a confidence level.
#'   Profiles will be computed up to a deviance
#'   of \code{\link{qchisq}(max_level, df = 1)}.
#' @param grid_len
#'   A positive integer. Step sizes chosen adaptively by
#'   \code{\link[TMB]{tmbprofile}} will generate approximately
#'   this many points on each side of a profile's minimum point.
#' @param parallel
#'   An \code{"\link{egf_parallel}"} object defining parallelization
#'   options.
#' @param trace
#'   A \link{logical} flag.
#'   If \code{TRUE}, then basic tracing messages indicating progress
#'   are printed. Depending on \code{fitted$control$trace}, these may
#'   be mixed with optimization output.
#' @param .subset,.append
#'   Index vectors to be used (if non-\code{\link{NULL}}) in place of
#'   the result of evaluating \code{subset} and \code{append}.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
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
#' \code{A}, \code{x = fitted$best[fitted$nonrandom]}, and \code{max_level}
#' are retained as \link{attributes}.
#'
#' @seealso \code{\link{confint.egf_profile}}, \code{plot.egf_profile}
#' @export
#' @importFrom Matrix sparseMatrix KhatriRao
#' @importMethodsFrom Matrix diag
#' @importFrom methods as is
#' @importFrom stats vcov
#' @importFrom TMB tmbprofile openmp
#' @import parallel
profile.egf <- function(fitted,
                        which = NULL,
                        A = NULL,
                        top = egf_get_names_top(fitted, link = TRUE),
                        subset = NULL,
                        append = NULL,
                        max_level = 0.99,
                        grid_len = 12,
                        parallel = egf_parallel(),
                        trace = TRUE,
                        .subset = NULL,
                        .append = NULL,
                        ...) {
  stopifnot(inherits(parallel, "egf_parallel"))
  stop_if_not_number_in_interval(max_level, 0, 1, "()")
  stop_if_not_number_in_interval(grid_len, 1, Inf, "[)")
  stop_if_not_true_false(trace)
  n <- length(fitted$nonrandom)

  ## If profiling user-specified elements of `c(beta, theta)`
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
  ## of elements of `c(beta, theta)`
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
    N <- sum(subset)
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
    r <- which
    h <- sqrt(diag(V)[which]) / 4
  } else {
    r <- lapply(seq_len(nrow(A)), function(i) A[i, ])
    ## Covariance matrix of `A %*% c(beta, theta)`
    V <- A %*% unclass(V) %*% t(A)
    h <- sqrt(diag(V)) / 4
  }
  ytol <- qchisq(max_level, df = 1) / 2 # y := diff(nll) = deviance / 2
  ystep <- ytol / grid_len

  tmbprofile_args <- list(
    obj = fitted$tmb_out,
    ytol = ytol,
    ystep = ystep,
    trace = FALSE
  )
  omp_num_threads <- fitted$control$omp_num_threads

  do_profile <- function(r, h, i) {
    if (trace) {
      cat(sprintf("Computing likelihood profile %d of %d...\n", i, m))
    }
    on <- openmp(n = NULL)
    if (on > 0L) {
      openmp(n = omp_num_threads)
      on.exit(openmp(n = on))
    }
    tmbprofile_args[[switch(method, which = "name", "lincomb")]] <- r
    tmbprofile_args$h <- h
    res <- do.call(tmbprofile, tmbprofile_args)
    i_min <- which.min(res[[2L]])
    res[[2L]] <- 2 * (res[[2L]] - res[i_min, 2L]) # deviance = 2 * diff(nll)
    names(res) <- c("value", "deviance")
    res[-i_min, , drop = FALSE] # `tmbprofile` duplicates this row
  }

  if (parallel$method == "snow") {
    environment(do_profile) <- .GlobalEnv # see comment in R/boot.R
    if (is.null(parallel$cl)) {
      cl <- do.call(makePSOCKcluster, parallel$options)
      on.exit(stopCluster(cl))
    } else {
      cl <- parallel$cl
    }
    clusterEvalQ(cl, library("TMB"))
    clusterExport(cl,
      varlist = c("tmbprofile_args", "trace", "m", "method", "omp_num_threads"),
      envir = environment()
    )
    res <- clusterMap(cl, do_profile, r = r, h = h, i = seq_len(m))
  } else {
    if (nzchar(parallel$outfile)) {
      outfile <- file(parallel$outfile, open = "wt")
      sink(outfile, type = "output")
      sink(outfile, type = "message")
    } else {
      res <- switch(parallel$method,
        multicore = do.call(mcMap, c(list(f = do_profile, r = r, h = h, i = seq_len(m)), parallel$options)),
        serial = Map(do_profile, r = r, h = h, i = seq_len(m))
      )
    }
    if (nzchar(parallel$outfile)) {
      sink(type = "output")
      sink(type = "message")
    }
  }

  nrow_res <- vapply(res, nrow, 0L)
  res <- data.frame(
    linear_combination = rep.int(gl(m, 1L), nrow_res),
    do.call(rbind, res),
    row.names = NULL
  )
  if (method == "top") {
    i <- rep.int(rep.int(subset, p), nrow_res)
    res <- data.frame(
      top = rep.int(rep.int(factor(top, levels = names_top), N), nrow_res),
      fitted$frame_windows[i, c("ts", "window"), drop = FALSE],
      res,
      combined[i, append, drop = FALSE],
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  attr(out, "A") <- A
  attr(out, "x") <- fitted$best[fitted$nonrandom]
  attr(out, "max_level") <- max_level
  attr(out, "method") <- method
  class(out) <- c("egf_profile", "data.frame")
  out
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
#' (level of \code{object$linear_combination}),
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
#' @export
#' @importFrom stats qchisq approx
confint.egf_profile <- function(object, parm, level = 0.95, link = TRUE, ...) {
  stop_if_not_true_false(link)
  stop_if_not_number_in_interval(level, 0, 1, "()")
  stopifnot(attr(object, "max_level") > level)
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
#' @param ...
#'   Optional graphical parameters passed to \code{\link{plot}},
#'   such as \code{type = "o"}. Note that \code{axes = FALSE} and
#'   \code{ann = FALSE} are hard-coded, so axes and axis titles
#'   may not be modifiable.
#'
#' @details
#' See topic \code{\link{egf_eval}} for details on nonstandard evaluation
#' of \code{subset}.
#'
#' @return
#' \code{\link{NULL}} (invisibly).
#'
#' @export
#' @import graphics
#' @importFrom stats confint
plot.egf_profile <- function(x, subset = NULL, sqrt = FALSE, level = NULL, ...) {
  subset <- egf_eval_subset(substitute(subset), x, parent.frame())
  m <- match(levels(factor(x$linear_combination[subset])), levels(x$linear_combination))

  stop_if_not_true_false(sqrt)
  f <- if (sqrt) base::sqrt else identity
  ymax <- f(max(x$deviance, na.rm = TRUE))
  ylab <- if (sqrt) expression(sqrt("deviance")) else "deviance"

  method <- attr(x, "method")
  do_ann_with_top <- (method == "top")

  do_segments <-
    !sqrt &&
    is.numeric(level) &&
    length(level) > 0L &&
    any(ok <- !is.na(level) & level > 0 & level < 1)
  if (do_segments) {
    level <- level[ok]
    ## Line segments at heights `h` in all plots
    h <- qchisq(level, df = 1)
    ## Line segment `j` to start at `v_lower[[i]][j]`
    ## and end at `v_upper[[i]][j]` in plot `i`
    ci <- lapply(level, function(p) confint(x, level = p))
    v_lower <- lapply(m, function(i) vapply(ci, `[`, 0, i, "lower"))
    v_upper <- lapply(m, function(i) vapply(ci, `[`, 0, i, "upper"))
  }

  op <- par(
    mar = c(3.5, 4, 1, 1),
    tcl = -0.4,
    cex.axis = 0.8,
    cex.lab = 0.9
  )
  on.exit(par(op))

  for (i in m) {
    r <- which(as.integer(x$linear_combination) == i)
    plot(
      f(deviance) ~ value,
      data = x,
      subset = r,
      ylim = c(0, ymax),
      ann = FALSE,
      axes = FALSE,
      ...
    )
    usr <- par("usr")
    if (do_segments) {
      v_lower <- vapply(ci, `[`, 0, i, "lower")
      v_upper <- vapply(ci, `[`, 0, i, "upper")
      segments(
        x0 = v_lower,
        x1 = v_upper,
        y0 = h,
        y1 = h,
        lty = 3
      )
      segments(
        x0 = c(v_lower, v_upper),
        x1 = c(v_lower, v_upper),
        y0 = usr[3L],
        y1 = rep.int(h, 2L),
        lty = 3
      )
      text(
        x = mean(usr[1:2]),
        y = h,
        labels = sprintf("%.3g%%", 100 * level),
        pos = 3,
        offset = 0.1,
        cex = 0.8
      )
    }
    box()
    axis(side = 1, mgp = c(3, 0.5, 0))
    axis(side = 2, mgp = c(3, 0.7, 0), las = 1)
    if (do_ann_with_top) {
      xlab <- as.character(x$top[r[1L]])
      main <- sprintf("window = %s", x$window[r[1L]])
      title(main = main, line = 2)
    } else {
      xlab <- sprintf("linear combination %d", i)
    }
    title(xlab = xlab, line = 2)
    title(ylab = ylab, line = 2.25)
  }

  invisible(NULL)
}
