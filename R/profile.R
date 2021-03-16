#' Compute likelihood profiles
#'
#' Computes the univariate likelihood profile of fixed effects
#' coefficients, log standard deviations of random effects
#' coefficients, and linear combinations thereof.
#'
#' @param fitted
#'   An `"egf"` object returned by [egf()].
#' @param which
#'   An integer vector indexing coefficients in `fitted$best`
#'   to be profiled. Must be a subset of `fitted$nonrandom`.
#' @param A
#'   A numeric matrix with `length(fitted$nonrandom)` columns.
#'   Each row specifies a linear combination of the elements of
#'   `fitted$best[fitted$nonrandom]` to be profiled. Ignored if
#'   `which` is non-`NULL`.
#' @param par
#'   A subset of `get_par_names(fitted, link = TRUE)` naming nonlinear
#'   model parameters whose population fitted values should be profiled.
#'   Ignored if `which` or `A` is non-`NULL`.
#' @param subset
#'   An expression to be evaluated in the combined model frame. Must
#'   evaluate to a logical vector or list of logical vectors indexing
#'   rows of the model frame, and thus fitting windows. Only population
#'   fitted values for the indexed fitting windows are profiled. The
#'   default (`NULL`) is to consider all fitting windows. Ignored if
#'   `which` or `A` is non-`NULL`.
#' @param append
#'   An expression indicating variables in the combined model frame
#'   to be included with the result. The default (`NULL`) is to
#'   append nothing. Ignored if `which` or `A` is non-`NULL`.
#' @param max_level
#'   A number in the interval (0,1). Profiles will be computed up to
#'   a deviance of `qchisq(max_level, df = 1)`.
#' @param grid_len
#'   A positive integer. Step sizes chosen adaptively by
#'   [TMB::tmbprofile()] will generate approximately this
#'   many points on each side of a profile's minimum point.
#' @param trace
#'   A logical scalar. If `TRUE`, then basic tracing messages
#'   are printed.
#' @inheritParams check_parallel
#' @inheritParams fitted.egf
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' The population fitted value is the fixed effects component of the
#' fitted value (see [fitted.egf()]). In fixed effects models, there
#' is no distinction.
#'
#' `which` is mapped to an `A` matrix composed of unit row vectors,
#' with 1 at array index `[i, which[i]]` for all `i`.
#' `A %*% fitted$best[fitted$nonrandom]` is precisely
#' `fitted$best[which]`.
#'
#' `par` and `subset` are mapped to a block `A` matrix. Block `k`
#' (from the top) is composed of rows of the fixed effects design
#' matrix for parameter `par[k]` (those rows indexed by `subset`).
#' `A %*% fitted$par[fitted$nonrandom]` is a vector listing the
#' population fitted values for each parameter named in `par`,
#' for each fitting window indexed by `subset`.
#'
#' @inheritSection fitted.egf Nonstandard evaluation
#' @inheritSection fitted.egf Warning
#'
#' @return
#' A data frame inheriting from class `"egf_profile"`, with variables:
#' \item{`linear_combination`}{
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
#'   (`par`-based calls only.)
#'   Nonlinear model parameter,
#'   from `get_par_names(fitted, link = TRUE)`.
#' }
#' \item{`ts`}{
#'   (`par`-based calls only.)
#'   Time series, from `levels(fitted$frame_ts$ts)`.
#' }
#' \item{`window`}{
#'   (`par`-based calls only.)
#'   Fitting window, from `levels(fitted$frame_ts$window)`.
#' }
#' `A` and `x = fitted$best[fitted$nonrandom]`
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
                        par = get_par_names(fitted, link = TRUE),
                        subset = NULL,
                        append = NULL,
                        max_level = 0.99,
                        grid_len = 12,
                        trace = TRUE,
                        parallel = c("serial", "multicore", "snow"),
                        cores = getOption("egf.cores", 2L),
                        outfile = NULL,
                        cl = NULL,
                        .subset = NULL,
                        .append = NULL,
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
    if (is.numeric(A) && is.null(dim(A))) {
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
  } else if (!is.null(par)) {
    method <- "par"
    frame <- do.call(cbind, fitted$frame_par)
    frame <- frame[!duplicated(names(frame))]

    par <- unique(match.arg(par, several.ok = TRUE))
    subset <- subset_to_index(substitute(subset), frame, parent.frame(),
                              .subset = .subset)
    append <- append_to_index(substitute(append), frame, parent.frame(),
                              .append = .append)

    p <- length(par)
    w <- length(subset)
    m <- p * w
    A <- rep.int(0, m * n)
    dim(A) <- c(m, n)
    for (k in seq_len(p)) {
      i <- (k - 1L) * w + seq_len(w)
      j <- which(fitted$tmb_args$data$X_info$par == par[k])
      A[i, j] <- fitted$tmb_args$data$X[subset, j]
    }

  ## Otherwise
  } else {
    stop("One of `A`, `which`, and `par` must be non-NULL.")
  }

  if (method == "which") {
    rl <- which
    ## Covariance matrix of `c(beta, log_sd_b)`
    vc <- vcov(fitted, full = TRUE)
    hl <- sqrt(diag(vc)[which]) / 4
  } else {
    rl <- lapply(seq_len(nrow(A)), function(i) A[i, ])
    ## Covariance matrix of `A %*% c(beta, log_sd_b)`
    vc <- A %*% vcov(fitted, full = TRUE) %*% t(A)
    hl <- sqrt(diag(vc)) / 4
  }
  ytol <- qchisq(max_level, df = 1) / 2 # y := diff(nll) = deviance / 2
  ystep <- ytol / grid_len

  a <- list(
    obj = fitted$tmb_out,
    ytol = ytol,
    ystep = ystep,
    trace = FALSE
  )
  i_of_m <- function(i, m = m) {
    sprintf("%*d of %d", nchar(m), i, m)
  }

  do_profile <- function(r, h, i) {
    if (trace) {
      cat("Computing likelihood profile", i_of_m(i), "...\n")
    }
    a[[switch(method, which = "name", "lincomb")]] <- r
    a$h <- h
    d <- do.call(tmbprofile, a)
    i_min <- which.min(d[[2L]])
    d[[2L]] <- 2 * (d[[2L]] - d[i_min, 2L]) # deviance = 2 * diff(nll)
    names(d) <- c("value", "deviance")
    d[-i_min, , drop = FALSE] # `tmbprofile()` duplicates this row
  }

  if (parallel == "snow") {
    environment(do_profile) <- environment(i_of_m) <- .GlobalEnv # see R/boot.R
    if (is.null(cl)) {
      if (is.null(outfile)) {
        outfile <- ""
      }
      cl <- makePSOCKcluster(cores, outfile = outfile)
      on.exit(stopCluster(cl))
    }
    clusterEvalQ(cl, library("TMB"))
    clusterExport(cl, varlist = c("a", "i_of_m", "m", "method"), envir = environment())
    dl <- clusterMap(cl, do_profile, r = rl, h = hl, i = seq_len(m))
  } else {
    if (!is.null(outfile)) {
      sink(outfile, type = "output")
      sink(outfile, type = "message")
    }
    if (parallel == "muticore") {
      dl <- mcmapply(do_profile, r = rl, h = hl, i = seq_len(m),
                     SIMPLIFY = FALSE, mc.cores = cores)
    } else { # "serial"
      dl <- Map(do_profile, r = rl, h = hl, i = seq_len(m))
    }
    if (!is.null(outfile)) {
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
    ts <- fitted$frame_ts$ts
    window <- fitted$frame_ts$window
    k <- !is.na(window) & !duplicated(window)
    pn <- get_par_names(fitted, link = TRUE)

    out <- data.frame(
      out,
      par = rep.int(rep.int(factor(par, levels = pn), w), dl_nrow),
      ts = rep.int(rep.int(ts[k][subset], p), dl_nrow),
      window = rep.int(rep.int(window[k][subset], p), dl_nrow),
      frame[rep.int(rep.int(subset, p), dl_nrow), append, drop = FALSE],
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  structure(out,
    class = c("egf_profile", "data.frame"),
    A = A,
    x = fitted$best[fitted$nonrandom],
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
#'   A number in the interval (0,1). The desired confidence level.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' Each supplied likelihood profile
#' (level of `object$linear_combination`),
#' is linearly interpolated to approximate the two solutions of
#' `deviance(value) = qchisq(level, df = 1)`.
#' These provide the lower and upper confidence limits of interest
#' (see Wilks' theorem).
#'
#' @return
#' A data frame with variables:
#' \item{`linear_combination`}{
#'   Row index of linear combination,
#'   from `seq_len(nrow(attr(object, "A")))`.
#' }
#' \item{`estimate`, `lower`, `upper`}{
#'   Estimate of linear combination
#'   and approximate lower and upper confidence limits.
#' }
#' `level`, `A = attr(object, "A")`, and `x = attr(object, "x")`
#' are retained as attributes.
#'
#' @export
#' @importFrom stats qchisq approx
confint.egf_profile <- function(object, parm, level = 0.95, ...) {
  stop_if_not_number_in_interval(level, 0, 1, "()")

  q <- qchisq(level, df = 1)
  stop_if_not(
    tapply(object$deviance, object$linear_combination, max, na.rm = TRUE) > q,
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
      d[1L, -(1:3), drop = FALSE],
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, lapply(split(object, object$linear_combination), f))
  out$linear_combination <- as.integer(as.character(out$linear_combination))
  row.names(out) <- NULL
  structure(out,
    level = level,
    A = attr(object, "A"),
    x = attr(object, "x")
  )
}

#' Plot likelihood profiles
#'
#' A method for inspecting computed likelihood profiles.
#'
#' @param x
#'   An `"egf_profile"` object returned by [profile.egf()].
#' @param subset
#'   An expression to be evaluated in `x`. Must evaluate to a
#'   logical vector or list of logical vectors indexing rows of
#'   `x`. Only indexed profiles are plotted. The default (`NULL`)
#'   is to plot all profiles.
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
  subset <- subset_to_index(substitute(subset), x, parent.frame())
  g <- factor(x$linear_combination,
    levels = levels(droplevels(x$linear_combination[subset]))
  )
  x_split <- split(x, g)

  stop_if_not_true_false(sqrt)
  f <- if (sqrt) base::sqrt else identity
  ymax <- f(max(x$deviance, na.rm = TRUE))
  ylab <- if (sqrt) expression(sqrt("deviance")) else "deviance"
  ann_with_par <- length(x) > 3L

  any_segments <- !sqrt && is.numeric(level) && length(level) > 0L
  if (any_segments) {
    stop_if_not(
      level > 0,
      level < 1,
      m = "Elements of `level` must be numbers\nin the interval (0,1)."
    )

    ## Line segment `j` at height `h[j]` in all plots
    h <- f(qchisq(level, df = 1))
    ## Line segment `j` to start at `v_lower[[i]][j]`
    ## and end at `v_upper[[i]][j]` in plot `i`
    cil <- lapply(level, function(p) confint(x, level = p))
    m <- match(levels(g), levels(x$linear_combination))
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
    if (any_segments) {
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
      xlab <- sprintf("linear combination %d", x_split[[i]]$linear_combination[1L])
    }
    title(xlab = xlab, line = 2)
    title(ylab = ylab, line = 2.25)
  }

  invisible(NULL)
}
