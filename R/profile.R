#' Compute likelihood profiles
#'
#' Computes the univariate likelihood profile of fixed effects
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
#' @param par
#'   A subset of \code{\link{get_par_names}(fitted, link = TRUE)}
#'   naming nonlinear and dispersion model parameters whose population
#'   fitted values (see \code{\link{fitted.egf}}) should be profiled.
#'   Ignored if \code{which} or \code{A} is non-\code{\link{NULL}}.
#' @param subset
#'   An expression to be evaluated in the combined model frame
#'   (see \code{\link{make_combined}}). Must evaluate to
#'   a \link{logical} vector indexing rows of the data frame,
#'   and thus fitting windows. Only population fitted values for
#'   indexed windows are profiled. The default (\code{\link{NULL}})
#'   is to consider all windows. Ignored if \code{which} or \code{A}
#'   is non-\code{\link{NULL}}.
#' @param append
#'   An expression indicating variables in the combined model frame
#'   (see \code{\link{make_combined}}) to be included with the result.
#'   The default (\code{\link{NULL}}) is to append nothing.
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
#'   An \code{"\link{egf_parallel}"} object defining parallelization options.
#' @param trace
#'   A \link{logical} flag. If \code{TRUE}, then basic tracing messages
#'   indicating progress are printed. Depending on \code{fitted$control$trace},
#'   these may be mixed with optimization output.
#' @param .subset
#'   A \link{logical} vector to be used (if non-\code{\link{NULL}})
#'   in place of the result of evaluating \code{subset}.
#' @param .append
#'   A \link{character} vector listing variable names to be used
#'   (if non-\code{\link{NULL}}) in place of the result of evaluating
#'   \code{append}.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' \code{which} is mapped to an \code{A} matrix composed of unit row vectors,
#' with 1 at array index \code{[i, which[i]]} for all \code{i}, such that
#' \code{A \link[=matmult]{\%*\%} fitted$best[fitted$nonrandom]} is precisely
#' \code{fitted$best[which]}.
#'
#' \code{par} and \code{subset} are mapped to a block \code{A} matrix.
#' Block \code{k} (from the top) is composed of rows of the fixed effects
#' design matrix for parameter \code{par[k]} (those rows indexed by
#' \code{subset}). \code{A \link[=matmult]{\%*\%} fitted$par[fitted$nonrandom]}
#' is a vector listing the population fitted values for each parameter
#' named in \code{par}, for each fitting window indexed by \code{subset}.
#'
#' See topic \code{\link{nse}} for details on nonstandard evaluation
#' of \code{subset} and \code{append}.
#'
#' @return
#' A \link[=data.frame]{data frame} inheriting from \link{class}
#' \code{"egf_profile"}, with variables:
#' \item{par}{
#'   (\code{par}-based calls only.)
#'   Nonlinear or dispersion model parameter,
#'   from \code{\link{get_par_names}(fitted, link = link)}.
#' }
#' \item{ts}{
#'   (\code{par}-based calls only.)
#'   Time series, from \code{\link{levels}(fitted$endpoints$ts)}.
#' }
#' \item{window}{
#'   (\code{par}-based calls only.)
#'   Fitting window, from \code{\link{levels}(fitted$endpoints$window)}.
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
#' @importFrom TMB tmbprofile
#' @import parallel
profile.egf <- function(fitted,
                        which = NULL,
                        A = NULL,
                        par = get_par_names(fitted, link = TRUE),
                        subset = NULL,
                        append = NULL,
                        max_level = 0.99,
                        grid_len = 12,
                        parallel = egf_parallel(),
                        trace = TRUE,
                        .subset = NULL,
                        .append = NULL,
                        ...) {
  stop_if_not(
    inherits(parallel, "egf_parallel"),
    m = "`parallel` must inherit from class \"egf_parallel\". See `?egf_parallel`."
  )
  stop_if_not_number_in_interval(max_level, 0, 1, "()")
  stop_if_not_number_in_interval(grid_len, 1, Inf, "[)")
  stop_if_not_true_false(trace)
  n <- length(fitted$nonrandom)

  ## If profiling user-specified elements of `c(beta, theta)`
  if (!is.null(which)) {
    method <- "which"
    stop_if_not(
      is.numeric(which),
      (which <- unique(which)) %in% fitted$nonrandom,
      m = "`which` must be a subset of `fitted$nonrandom`."
    )
    m <- length(which)
    A <- sparseMatrix(i = seq_len(m), j = which, x = 1, dims = c(m, n))

  ## If profiling user-specified linear combinations
  ## of elements of `c(beta, theta)`
  } else if (!is.null(A)) {
    method <- "A"
    if (is.numeric(A) && is.null(dim(A))) {
      dim(A) <- c(1L, length(A))
    }
    stop_if_not(
      (is.matrix(A) && is.numeric(A)) || is(A, "dMatrix"),
      all(is.finite(A)),
      m = "`A` must be a finite, numeric matrix or vector."
    )
    stop_if_not(
      ncol(A) == n,
      m = "`A` must have `length(fitted$nonrandom)` columns."
    )
    stop_if_not(
      rowSums(abs(A)) > 0,
      m = "`A` must have at least one nonzero element in each row."
    )
    m <- nrow(A)

  ## If profiling population fitted values
  ## of nonlinear and dispersions model parameters
  } else if (!is.null(par)) {
    method <- "par"
    par_names <- get_par_names(fitted, link = TRUE)
    par <- unique(match.arg(par, par_names, several.ok = TRUE))
    combined <- make_combined(fitted)
    subset <- subset_to_index(substitute(subset), data = combined, enclos = parent.frame(),
                              .subset = .subset)
    append <- append_to_index(substitute(append), data = combined, enclos = parent.frame(),
                              .append = .append)

    p <- length(par)
    N <- sum(subset)
    m <- p * N

    f <- factor(fitted$tmb_args$data$X_info$par, levels = par)
    J <- as(f, "sparseMatrix")
    X <- fitted$tmb_args$data$X[subset, , drop = FALSE]
    A <- KhatriRao(J, X)
    if (has_random(fitted)) {
      index <- grepl("^beta\\[", names(fitted$best)[fitted$nonrandom])
      A@p <- c(0L, replace(rep_len(NA_integer_, n), index, A@p[-1L]))
      A@p <- locf(A@p)
      A@Dim <- c(m, n)
    }

  ## Otherwise
  } else {
    stop("One of `A`, `which`, and `par` must be non-NULL.")
  }

  ## Covariance matrix of `c(beta, theta)`
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

  do_profile <- function(r, h, i) {
    if (trace) {
      cat(sprintf("Computing likelihood profile %*d of %d...\n", nchar(m), i, m))
    }
    tmbprofile_args[[switch(method, which = "name", "lincomb")]] <- r
    tmbprofile_args$h <- h
    d <- do.call(tmbprofile, tmbprofile_args)
    i_min <- which.min(d[[2L]])
    d[[2L]] <- 2 * (d[[2L]] - d[i_min, 2L]) # deviance = 2 * diff(nll)
    names(d) <- c("value", "deviance")
    d[-i_min, , drop = FALSE] # `tmbprofile` duplicates this row
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
    clusterExport(cl, varlist = c("tmbprofile_args", "m", "method"), envir = environment())
    dl <- clusterMap(cl, do_profile, r = r, h = h, i = seq_len(m))

  } else {
    if (parallel$outfile != "") {
      outfile <- file(parallel$outfile, open = "wt")
      sink(outfile, type = "output")
      sink(outfile, type = "message")
    } else {
      dl <- switch(parallel$method,
        multicore = do.call(mcMap, c(list(f = do_profile, r = r, h = h, i = seq_len(m)), parallel$options)),
        serial = Map(do_profile, r = r, h = h, i = seq_len(m))
      )
    }
    if (parallel$outfile != "") {
      sink(type = "output")
      sink(type = "message")
    }
  }

  dl_nrow <- vapply(dl, nrow, 0L)
  out <- data.frame(
    linear_combination = rep.int(gl(m, 1L), dl_nrow),
    do.call(rbind, dl),
    row.names = NULL
  )
  if (method == "par") {
    w <- rep.int(rep.int(which(subset), p), dl_nrow)
    out <- data.frame(
      par = rep.int(rep.int(factor(par, levels = par_names), N), dl_nrow),
      fitted$endpoints[w, c("ts", "window"), drop = FALSE],
      out,
      combined[w, append, drop = FALSE],
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
#'   of nonlinear or dispersion model parameters,
#'   then confidence intervals on inverse link-transformed
#'   fitted values are returned.
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
#'   if \code{link = FALSE} and linear combinations represent
#'   population fitted values of nonlinear or dispersion model
#'   parameters.
#' }
#' \link[=attributes]{Attributes} \code{A} and \code{x} of \code{object}
#' are retained as attributes of the result. So is \code{level}.
#'
#' @export
#' @importFrom stats qchisq approx
confint.egf_profile <- function(object, parm, level = 0.95, link = TRUE, ...) {
  stop_if_not_true_false(link)
  stop_if_not_number_in_interval(level, 0, 1, "()")
  stop_if_not(
    attr(object, "max_level") > level,
    m = "`attr(object, \"max_level\")` must exceed `level`."
  )

  s <- c("par", "ts", "window", "linear_combination", "value", "deviance")
  if (attr(object, "method") == "par") {
    k1 <- 1:4
    k2 <- -match(s, names(object), 0L)
  } else {
    k1 <- 4L
    k2 <- -match(s[4:6], names(object), 0L)
  }
  q <- qchisq(level, df = 1)

  do_solve <- function(d) {
    i_min <- which.min(d$deviance)
    i_left <- seq_len(i_min)
    i_right <- seq.int(i_min, nrow(d))
    estimate <- d$value[i_min]
    lower <- approx(x = d$deviance[i_left],  y = d$value[i_left],  xout = q)$y
    upper <- approx(x = d$deviance[i_right], y = d$value[i_right], xout = q)$y
    data.frame(
      d[1L, k1, drop = FALSE],
      estimate,
      lower,
      upper,
      d[1L, k2, drop = FALSE],
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, by(object, object$linear_combination, do_solve, simplify = FALSE))
  if (attr(object, "method") == "par" && !link) {
    elu <- c("estimate", "lower", "upper")
    out[elu] <- lpapply(out[elu], out$par,
      f = lapply(string_extract_link(levels(out$par)), match_link, inverse = TRUE)
    )
    levels(out$par) <- string_remove_link(levels(out$par))
  }
  out$linear_combination <- as.integer(as.character(out$linear_combination))
  row.names(out) <- NULL
  attr(out, "A") <- attr(object, "A")
  attr(out, "x") <- attr(object, "x")
  attr(out, "level") <- level
  out
}

# #' Plot likelihood profiles
# #'
# #' A method for inspecting computed likelihood profiles.
# #'
# #' @param x
# #'   An `"egf_profile"` object returned by [profile.egf()].
# #' @param subset
# #'   An expression to be evaluated in `x`. Must evaluate to a
# #'   logical vector indexing rows of `x`. Only indexed profiles
# #'   are plotted. The default (`NULL`) is to plot all profiles.
# #' @param sqrt
# #'   A logical scalar. If `TRUE`, then square root-transformed
# #'   deviance is plotted.
# #' @param level
# #'   A numeric vector with elements in (0,1), or otherwise `NULL`.
# #'   If non-`NULL` and `sqrt = FALSE`, then line segments are
# #'   drawn showing the intersection of the profile with lines
# #'   `deviance = qchisq(level, df = 1)`.
# #' @param ...
# #'   Optional graphical parameters passed to [graphics::plot()],
# #'   such as `type = "o"`. Note that `axes = FALSE` and `ann = FALSE`
# #'   are hard-coded, so axes and axis titles cannot be modified.
# #'
# #' @details
# #' See topic [`nse`] for details on nonstandard evaluation of `subset`.
# #'
# #' @return
# #' `NULL` (invisibly).
# #'
# #' @export
# #' @import graphics
# #' @importFrom stats confint
# plot.egf_profile <- function(x, subset = NULL, sqrt = FALSE,
#                              level = NULL, ...) {
#   subset <- subset_to_index(substitute(subset), x, parent.frame())
#   g <- factor(x$linear_combination,
#     levels = levels(droplevels(x$linear_combination[subset]))
#   )
#   x_split <- split(x, g)
#
#   stop_if_not_true_false(sqrt)
#   f <- if (sqrt) base::sqrt else identity
#   ymax <- f(max(x$deviance, na.rm = TRUE))
#   ylab <- if (sqrt) expression(sqrt("deviance")) else "deviance"
#   ann_with_par <- length(x) > 3L
#
#   any_segments <- !sqrt && is.numeric(level) && length(level) > 0L
#   if (any_segments) {
#     stop_if_not(
#       level > 0,
#       level < 1,
#       m = "Elements of `level` must be numbers\nin the interval (0,1)."
#     )
#
#     ## Line segment `j` at height `h[j]` in all plots
#     h <- f(qchisq(level, df = 1))
#     ## Line segment `j` to start at `v_lower[[i]][j]`
#     ## and end at `v_upper[[i]][j]` in plot `i`
#     cil <- lapply(level, function(p) confint(x, level = p))
#     m <- match(levels(g), levels(x$linear_combination))
#     v_lower <- lapply(m, function(i) vapply(cil, `[`, 0, i, "lower"))
#     v_upper <- lapply(m, function(i) vapply(cil, `[`, 0, i, "upper"))
#   }
#
#   op <- par(
#     mar = c(3.5, 4, 1, 1),
#     tcl = -0.4,
#     cex.axis = 0.8,
#     cex.lab = 0.9
#   )
#   on.exit(par(op))
#
#   for (i in seq_along(x_split)) {
#     plot(
#       f(deviance) ~ value,
#       data = x_split[[i]],
#       ylim = c(0, ymax),
#       ann = FALSE,
#       axes = FALSE,
#       ...
#     )
#     if (any_segments) {
#       segments(
#         x0 = v_lower[[i]],
#         x1 = v_upper[[i]],
#         y0 = h,
#         y1 = h,
#         lty = 3
#       )
#       segments(
#         x0 = c(v_lower[[i]], v_upper[[i]]),
#         x1 = c(v_lower[[i]], v_upper[[i]]),
#         y0 = par("usr")[3L],
#         y1 = rep.int(h, 2L),
#         lty = 3
#       )
#       text(
#         x = mean(par("usr")[1:2]),
#         y = h,
#         labels = sprintf("%.3g%%", 100 * level),
#         pos = 3, offset = 0.1, cex = 0.8
#       )
#     }
#     box()
#     axis(side = 1, mgp = c(3, 0.5, 0))
#     axis(side = 2, mgp = c(3, 0.7, 0), las = 1)
#     if (ann_with_par) {
#       xlab <- as.character(x_split[[i]]$par[1L])
#       main <- sprintf("window = %s", x_split[[i]]$window[1L])
#       title(main = main, line = 2)
#     } else {
#       xlab <- sprintf("linear combination %d", x_split[[i]]$linear_combination[1L])
#     }
#     title(xlab = xlab, line = 2)
#     title(ylab = ylab, line = 2.25)
#   }
#
#   invisible(NULL)
# }
