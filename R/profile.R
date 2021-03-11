#' Compute likelihood profiles
#'
#' Computes the univariate likelihood profile of fixed effects
#' coefficients, log standard deviations of random effects
#' coefficients, and linear combinations thereof.
#'
#' @param fitted
#'   An `"egf"` object returned by [egf()].
#' @param which
#'   An integer vector indexing coefficients in `fitted$par`
#'   to be profiled. Must be a subset of `fitted$nonrandom`.
#' @param A
#'   A numeric matrix with `length(fitted$nonrandom)` columns.
#'   Each row specifies a linear combination of the elements
#'   of `fitted$par[fitted$nonrandom]` to be profiled. Ignored
#'   if `which` is non-`NULL`.
#' @param parm
#'   A character vector listing nonlinear model parameters whose
#'   population fitted values should be profiled. Must be a subset
#'   of `get_par_names(fitted, link = TRUE)`. Ignored if one of
#'   `which` and `lin_comb` is non-`NULL`. (Here, "population
#'   fitted value" refers to the fixed effects component of the
#'   individual or overall fitted value, which accounts for both
#'   fixed and random effects.)
#' @param max_level
#'   A number in the interval (0,1). Profiles will be computed
#'   up to a deviance of `qchisq(max_level, df = 1)`.
#' @param grid_len
#'   A positive integer. Step sizes chosen adaptively by
#'   [TMB::tmbprofile()] will generate approximately this
#'   many points on each side of a profile's minimum point.
#' @param trace
#'   A logical scalar. If `TRUE`, then basic tracing messages
#'   are printed.
#' @inheritParams check_parallel
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' `which` is mapped to an `A` matrix composed of unit row vectors,
#' with 1 at array index `[i, which[i]]` for all `i`.
#' `A %*% fitted$par[fitted$nonrandom]` reproduces `fitted$par[which]`.
#'
#' `parm` is mapped to an `A` matrix composed of blocks that
#' are precisely the fixed effects design matrices for nonlinear
#' model parameters `parm`. `A %*% fitted$par[fitted$nonrandom]`
#' is a vector listing each parameters' population fitted values
#' (one per fitting window per parameter).
#'
#' @return
#' A data frame inheriting from class `"egf_profile"`, with variables:
#' \item{`index`}{
#'   Row index of linear combination, from `seq_len(nrow(A))`.
#' }
#' \item{`value`}{
#'   Value of linear combination being profiled.
#' }
#' \item{`deviance`}{
#'   Deviance of the restricted model that assumes `value`
#'   for the linear combination being profiled.
#' }
#' \item{`par`}{
#'   (`parm`-based calls only.)
#'   Nonlinear model parameter, from `parm`.
#' }
#' \item{`window`}{
#'   (`parm`-based calls only.)
#'   Fitting window, from `levels(fitted$frame_ts$.window)`.
#' }
#' Matrix `A` and vector `par = fitted$par[fitted$nonrandom]`
#' are retained as attributes.
#'
#' @seealso [confint.egf_profile()], [plot.egf_profile()]
#' @export
#' @importFrom stats vcov
#' @importFrom TMB tmbprofile
#' @import parallel
profile.egf <- function(fitted,
                        which = NULL,
                        A = NULL,
                        parm = get_par_names(fitted, link = TRUE),
                        max_level = 0.99,
                        grid_len = 12L,
                        trace = TRUE,
                        parallel = c("serial", "multicore", "snow"),
                        cores = getOption("egf.cores", 2L),
                        outfile = NULL,
                        cl = NULL,
                        ...) {
  stop_if_not_number_in_interval(max_level, 0, 1, "()")
  stop_if_not_number_in_interval(grid_len, 1, Inf, "[)")
  stop_if_not_true_false(trace)
  parallel <- match.arg(parallel)
  check_parallel(parallel, cores, outfile, cl)
  n <- length(fitted$nonrandom)

  ## If profiling user-specified elements of `c(beta, log_sd_b)`
  if (!is.null(which)) {
    method <- "which"
    stop_if_not(
      is.numeric(which),
      length(which) > 0L,
      (which <- unique(which)) %in% fitted$nonrandom,
      m = "`which` must be a subset of `fitted$nonrandom`."
    )
    m <- length(which)
    A <- rep.int(0L, m * n)
    dim(A) <- c(m, n)
    A[cbind(seq_len(m), which)] <- 1L

  ## If profiling user-specified linear combinations
  ## of elements of `c(beta, log_sd_b)`
  } else if (!is.null(A)) {
    method <- "A"
    if (is.vector(A)) {
      dim(A) <- c(1L, length(A))
    }
    stop_if_not(
      is.matrix(A),
      is.numeric(A),
      is.finite(A),
      m = "`A` must be a finite, numeric matrix or vector."
    )
    stop_if_not(
      nrow(A) > 0L,
      ncol(A) == n,
      m = "`A` must have at least one row and exactly\n`length(fitted$nonrandom)` columns."
    )
    stop_if_not(
      rowSums(abs(A)) > 0,
      m = "`A` must have at least one nonzero element in each row."
    )
    m <- nrow(A)

  ## If profiling population fitted values of nonlinear model parameters
  } else if (!is.null(parm)) {
    method <- "parm"
    parm <- unique(match.arg(parm, several.ok = TRUE))
    p <- length(parm)
    w <- nlevels(fitted$frame_ts$.window)
    m <- p * w
    A <- rep.int(0, m * n)
    dim(A) <- c(m, n)
    for (k in seq_len(p)) {
      i <- (k - 1L) * w + seq_len(w)
      j <- (fitted$tmb_args$data$X_info$par == parm[k])
      A[i, j] <- fitted$tmb_args$data$X[, j]
    }

  ## Otherwise
  } else {
    stop("One of `A`, `which`, and `parm` must be non-NULL.")
  }

  if (method == "which") {
    xl <- which
    ## Covariance matrix of `c(beta, log_sd_b)`
    vc <- vcov(fitted, full = TRUE)
    hl <- sqrt(diag(vc)[which]) / 4
  } else {
    xl <- lapply(seq_len(nrow(A)), function(i) A[i, ])
    ## Covariance matrix of `A %*% c(beta, log_sd_b)`
    vc <- A %*% vcov(fitted, full = TRUE) %*% t(A)
    hl <- sqrt(diag(vc)) / 4
  }
  ytol <- qchisq(max_level, df = 1) / 2 # y := diff(nll) = deviance / 2
  ystep <- ytol / grid_len

  a <- list(
    obj   = fitted$tmb_out,
    ytol  = ytol,
    ystep = ystep,
    trace = FALSE
  )
  i_of_m <- function(i, m = m) {
    sprintf("%*d of %d", nchar(m), i, m)
  }

  do_profile <- function(x, h, i) {
    if (trace) {
      cat("Computing likelihood profile", i_of_m(i), "...\n")
    }
    a[[switch(method, which = "name", "lincomb")]] <- x
    a[["h"]] <- h
    d <- do.call(tmbprofile, a)
    d[[2L]] <- 2 * (d[[2L]] - min(d[[2L]], na.rm = TRUE)) # deviance = 2 * diff(nll)
    names(d) <- c("value", "deviance")
    unique(d) # omits instances of zero step size
  }

  if (parallel == "snow") {
    ## See comment in R/boot.R
    environment(do_profile) <- environment(i_of_m) <- .GlobalEnv
    if (is.null(cl)) {
      if (is.null(outfile)) {
        outfile <- ""
      }
      cl <- makePSOCKcluster(cores, outfile = outfile)
      on.exit(stopCluster(cl))
    }
    clusterEvalQ(cl, library("TMB"))
    clusterExport(cl, varlist = c("a", "method", "i_of_m", "m"), envir = environment())
    dl <- clusterMap(cl, do_profile, x = xl, h = hl, i = seq_len(m))
  } else {
    if (!is.null(outfile)) {
      sink(outfile, type = "output")
      sink(outfile, type = "message")
    }
    dl <- switch(parallel,
      multicore = mcmapply(do_profile, x = xl, h = hl, i = seq_len(m),
                           SIMPLIFY = FALSE, mc.cores = cores),
      serial    =      Map(do_profile, x = xl, h = hl, i = seq_len(m))
    )
    if (!is.null(outfile)) {
      sink(type = "output")
      sink(type = "message")
    }
  }

  dl_nrow <- vapply(dl, nrow, 0L)
  out <- data.frame(
    index = rep.int(gl(m, 1L), dl_nrow),
    do.call(rbind, dl),
    row.names = NULL
  )
  if (method == "parm") {
    l <- levels(fitted$frame_ts$.window)
    out$par <- rep.int(gl(p, w, labels = parm), dl_nrow)
    out$window = rep.int(rep.int(gl(w, 1L, labels = l), p), dl_nrow)
  }
  structure(out,
    class = c("egf_profile", "data.frame"),
    A = A,
    par = fitted$par[fitted$nonrandom],
  )
}

#' Confidence intervals from likelihood profiles
#'
#' Computes confidence intervals on fixed effects coefficients,
#' log standard deviations of random effects coefficients, and
#' linear combinations thereof from their univariate likelihood
#' profiles.
#'
#' @param object
#'   An `"egf_profile"` object returned by [profile.egf()].
#' @param parm
#'   Unused argument included for generic consistency.
#' @param level
#'   A number in the interval (0,1), indicating a confidence level.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' For each supplied likelihood profile (level of `object$index`),
#' [stats::approx()] is called to approximate (via linear interpolation)
#' the two solutions of `deviance(value) = qchisq(level, df = 1)`,
#' which provide the lower and upper confidence limits of interest.
#'
#' @return
#' A data frame with variables:
#' \item{`index`}{
#'   Row index of linear combination, from `levels(object$index)`.
#' }
#' \item{`estimate`, `lower`, `upper`}{
#'   Estimate and (approximate) lower and upper confidence limits
#'   of linear combination.
#' }
#' \item{`par`}{
#'   (If `object$par` is non-`NULL`.)
#'   Nonlinear model parameter, from `levels(object$par)`.
#' }
#' \item{`window`}{
#'   (If `object$window` is non-`NULL`.)
#'   Fitting window, from `levels(object$window)`.
#' }
#' `level` and attributes `A` and `par` of `object` are retained
#' as attributes.
#'
#' @export
#' @importFrom stats qchisq approx
confint.egf_profile <- function(object, parm, level = 0.95, ...) {
  q <- qchisq(level, df = 1)
  stop_if_not_number_in_interval(level, 0, 1, "()")
  stop_if_not(
    tapply(object$deviance, object$index, max, na.rm = TRUE) > q,
    m = paste0(
      "Maximum deviance must exceed `qchisq(level, df = 1)`.\n",
      "Reprofile with higher `max_level` or retry with lower `level`."
    )
  )

  f <- function(d) {
    i_min <- which.min(d$deviance)
    i_left <- seq_len(i_min)
    i_right <- seq.int(i_min, nrow(d))

    estimate <- d$value[i_min]
    lower <- approx(
      x = d$deviance[i_left],
      y = d$value[i_left],
      xout = q
    )$y
    upper <- approx(
      x = d$deviance[i_right],
      y = d$value[i_right],
      xout = q
    )$y
    data.frame(
      d[1L, 1L, drop = FALSE],
      estimate,
      lower,
      upper,
      d[1L, -(1:3), drop = FALSE]
    )
  }

  out <- do.call(rbind, lapply(split(object, object$index), f))
  out$index <- as.integer(as.character(out$index))
  row.names(out) <- NULL
  structure(out,
    A = attr(object, "A"),
    par = attr(object, "par"),
    level = level
  )
}

#' Plot likelihood profiles
#'
#' A method for inspecting computed likelihood profiles.
#'
#' @param x
#'   An `"egf_profile"` object returned by [profile.egf()].
#' @param subset
#'   A list of logical expressions to be evaluated in `x[-(2:3)]`
#'   (`x` without variables `value` and `deviance`). Only indexed
#'   profiles are plotted. Use the default (`NULL`) to plot all
#'   profiles.
#' @param sqrt
#'   A logical scalar. If `TRUE`, then square root-transformed
#'   deviance is plotted.
#' @param level
#'   A numeric vector with elements in (0,1), or otherwise `NULL`.
#'   If non-`NULL` and `sqrt = FALSE`, then line segments are
#'   drawn showing the intersection of the profile with lines
#'   `deviance = qchisq(level, df = 1)`.
#' @param ...
#'   Optional graphical parameters passed to [graphics::plot()],
#'   such as `type = "o"`. Note that `axes = FALSE` and `ann = FALSE`
#'   are hard-coded, so axes and axis titles cannot be modified.
#'
#' @return
#' `NULL` (invisibly).
#'
#' @export
#' @import graphics
#' @importFrom stats confint
plot.egf_profile <- function(x, subset = NULL, sqrt = FALSE,
                             level = NULL, ...) {
  if (is.null(subset)) {
    g <- x$index
  } else {
    e <- substitute(subset)
    l <- eval(e, envir = x[-(2:3)], enclos = parent.frame())
    stop_if_not(
      is.list(l),
      vapply(l, is.vector, FALSE, "logical"),
      lengths(l) == nrow(x),
      m = sprintf("`subset` must evaluate to a list\nof logical vectors of length %d.", nrow(x))
    )
    lr <- Reduce(`&`, l)
    g <- factor(x$index, levels = levels(droplevels(x$index[lr & !is.na(lr)])))
  }
  x_split <- split(x, g)

  stop_if_not_true_false(sqrt)
  f <- if (sqrt) base::sqrt else identity
  ymax <- f(max(x$deviance, na.rm = TRUE))
  ylab <- if (sqrt) expression(sqrt("deviance")) else "deviance"
  ann_with_par <- "par" %in% names(x)

  if (!sqrt & !is.null(level)) {
    stop_if_not(
      is.vector(level, "numeric"),
      length(level) > 0L,
      level > 0,
      level < 1,
      m = "`level` must be NULL or a numeric vector\nwith elements in (0,1)."
    )
    ## Line segment `j` at height `h[j]` in all plots
    h <- f(qchisq(level, df = 1))
    ## Line segment `j` to start at `v_lower[[i]][j]`
    ## and end at `v_upper[[i]][j]` in plot `i`
    cil <- lapply(level, function(p) confint(x, level = p))
    m <- match(levels(g), levels(x$index))
    v_lower <- lapply(m, function(i) vapply(cil, `[`, 0, i, "lower"))
    v_upper <- lapply(m, function(i) vapply(cil, `[`, 0, i, "upper"))
  }

  op <- par(
    mar = c(3.5, 4, 1, 1),
    tcl = -0.4,
    cex.axis = 0.8,
    cex.lab = 0.9
  )
  on.exit(par(op))

  for (i in seq_along(x_split)) {
    plot(
      f(deviance) ~ value,
      data = x_split[[i]],
      ylim = c(0, ymax),
      ann = FALSE,
      axes = FALSE,
      ...
    )
    if (!sqrt && !is.null(level)) {
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
        y0 = par("usr")[3L],
        y1 = rep.int(h, 2L),
        lty = 3
      )
      text(
        x = mean(par("usr")[1:2]),
        y = h,
        labels = sprintf("%.3g%%", 100 * level),
        pos = 3, offset = 0.1, cex = 0.8
      )
    }
    box()
    axis(side = 1, mgp = c(3, 0.5, 0))
    axis(side = 2, mgp = c(3, 0.7, 0), las = 1)
    if (ann_with_par) {
      xlab <- as.character(x_split[[i]]$par[1L])
      main <- sprintf("window = %s", x_split[[i]]$window[1L])
      title(main = main, line = 2)
    } else {
      xlab <- sprintf("linear combination %d", x_split[[i]]$index[1L])
    }
    title(xlab = xlab, line = 2)
    title(ylab = ylab, line = 2.25)
  }

  invisible(NULL)
}
