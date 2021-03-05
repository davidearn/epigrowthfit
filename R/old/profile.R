#' Compute likelihood profiles
#'
#' Computes the (univariate) likelihood profile of fixed effect
#' coefficients, log standard deviations of random effects, and
#' linear combinations thereof.
#'
#' @param fitted
#'   An `"egf"` object returned by [egf()].
#' @param index
#'   An integer vector indexing elements of `fitted$par` to be profiled.
#'   Must be a subset of `fitted$nonrandom`.
#' @param lin_comb
#'   A numeric matrix with `length(fitted$nonrandom)` columns,
#'   whose rows specify linear combinations of the elements of
#'   `fitted$par[fitted$nonrandom]`. Ignored if `index` is non-`NULL`.
#' @param parm
#'   A character vector giving a subset of the response variables
#'   listed in `get_par_names(fitted, link = TRUE)`. The elements
#'   of `fitted$par[nonrandom]` belonging to the mixed effects model
#'   for these response variables, or linear combinations thereof
#'   (see `decontrast`), will be profiled. Ignored if `index` or
#'   `lin_comb` are non-`NULL`.
#' @param decontrast
#'   A logical scalar. If `TRUE`, then those linear combinations of
#'   the fixed effects coefficients (zero-sum contrasts) that recover
#'   within-group means are profiled instead of the fixed effects
#'   coefficients themselves. Ignored if `index` or `lin_comb` are
#'   non-`NULL`.
#' @param max_level
#'   A number in the interval (0,1). Profiles will be computed up to
#'   a deviance of `qchisq(max_level, df = 1)`.
#' @param grid_len
#'   A positive integer. Step sizes chosen adaptively by
#'   [TMB::tmbprofile()] will generate approximately this
#'   many points on each side of the minimum point.
#' @inheritParams boot_par
#'
#' @return
#' A data frame with variables:
#' \item{`name`}{
#'   A factor. If `index` or `parm` are used, then the name of the
#'   parameter or linear combination being profiled, taken or formed
#'   from `names(fitted$par)`. If `lin_comb` is used, then the factor
#'   levels are `seq_len(nrow(lin_comb))` corresponding to the rows
#'   of `lin_comb`.
#' }
#' \item{`value`}{
#'   A numeric vector. The value of the parameter
#'   or linear combination being profiled.
#' }
#' \item{`deviance`}{
#'   A numeric vector. The deviance of the restricted
#'   model that assumes `value` for the parameter or
#'   linear combination being profiled.
#' }
#'
#' @seealso [confint.egf_profile()], [plot.egf_profile()]
#' @export
#' @importFrom stats vcov
#' @import parallel
profile.egf <- function(fitted, index = NULL, lin_comb = NULL,
                        parm = get_par_names(fitted),
                        decontrast = FALSE,
                        max_level = 0.99, grid_len = 12,
                        trace = TRUE,
                        parallel = c("serial", "multicore", "snow"),
                        cores = getOption("egf.cores", 2L),
                        outfile = NULL,
                        cl = NULL,
                        ...) {
  stop_if_not(
    is.numeric(max_level),
    length(max_level) == 1L,
    max_level > 0,
    max_level < 1,
    m = "`max_level` must be a number in the interval (0,1)."
  )
  stop_if_not(
    is.numeric(grid_len),
    length(grid_len) == 1L,
    grid_len >= 1,
    !is.infinite(grid_len),
    m = "`grid_len` must be 1 or greater."
  )
  stop_if_not_tf(trace)
  parallel <- match.arg(parallel)
  check_parallel(parallel, cores, outfile, cl)

  ## If profiling user-specified elements of `u = c(beta, log_sd_b)`
  if (!is.null(index)) {
    stop_if_not(
      is.numeric(index),
      length(index) > 0L,
      index %in% fitted$nonrandom,
      m = "`index` must be a subset of `fitted$nonrandom`."
    )
    index <- unique(index)
    lin_comb <- NULL
    profile_names <- names(fitted$par)[index]

  ## If profiling user-specified linear combinations of elements of `u`
  } else if (!is.null(lin_comb)) {
    if (is.vector(lin_comb)) {
      lin_comb <- matrix(lin_comb, nrow = 1L)
    }
    stop_if_not(
      is.matrix(lin_comb),
      is.numeric(lin_comb),
      is.finite(lin_comb),
      m = "`lin_comb` must be a finite, numeric matrix or vector."
    )
    stop_if_not(
      nrow(lin_comb) > 0L,
      ncol(lin_comb) == length(fitted$nonrandom),
      apply(lin_comb, 1L, function(x) any(x != 0)),
      m = paste0(
        "`lin_comb` must have `length(fitted$nonrandom)` columns\n",
        "and at least one nonzero element per row."
      )
    )
    profile_names <- seq_len(nrow(lin_comb))

  ## If profiling (linear combinations of) elements of `u`
  ## relevant to user-specified response variables
  } else if (!is.null(parm)) {
    pn <- get_par_names(fitted, link = TRUE)
    stop_if_not(
      is.character(parm),
      length(parm) > 0L,
      (parm <- add_link_string(unique(remove_link_string(parm)))) %in% pn,
      m = "`parm` must be a subset of `get_par_names(fitted)`."
    )
    stop_if_not_tf(decontrast)

    ## Want relevant elements of `v = c(decontrast(beta), log_sd_b)`
    if (decontrast) {
      ## Matrix of linear coefficients such that
      ## `lin_comb %*% u = v[relevant]`
      lin_comb <- make_lin_comb_for_parm(fitted, parm)
      profile_names <- rownames(lin_comb)
    ## Want relevant elements of `u`
    } else {
      index <- make_index_for_parm(fitted, parm)
      lin_comb <- NULL
      profile_names <- names(index)
    }

  } else {
    stop("One of `index`, `lin_comb`, and `parm` must be non-NULL.")
  }

  if (is.null(lin_comb)) {
    ## Covariance matrix of `u`
    vc <- vcov(fitted, full = TRUE)
    hl <- sqrt(diag(vc)[index]) / 4
    xl <- unname(index)
  } else {
    ## Covariance matrix of `lin_comb %*% u`
    vc <- lin_comb %*% vcov(fitted, full = TRUE) %*% t(lin_comb)
    hl <- sqrt(diag(vc)) / 4
    xl <- lapply(seq_len(nrow(lin_comb)), function(i) lin_comb[i, ])
  }
  ytol <- qchisq(max_level, df = 1) / 2 # y := diff(nll) = deviance / 2
  ystep <- ytol / grid_len

  n <- length(xl)
  iofn <- function(i) sprintf("%*d of %d", nchar(n), i, n)
  f <- function(x, h, i) {
    if (trace) cat("Computing likelihood profile", iofn(i), "...\n")
    if (length(x) > 1L) {
      d <- TMB::tmbprofile(obj = fitted$tmb_out, lincomb = x, h = h,
                           ytol = ytol, ystep = ystep, trace = FALSE)
    } else {
      d <- TMB::tmbprofile(obj = fitted$tmb_out, name = x, h = h,
                           ytol = ytol, ystep = ystep, trace = FALSE)
    }
    d[[2L]] <- 2 * (d[[2L]] - min(d[[2L]], na.rm = TRUE)) # deviance = 2 * diff(nll)
    names(d) <- c("value", "deviance")
    d[!duplicated(d), , drop = FALSE] # omits instances of zero step size
  }

  if (parallel == "snow") {
    ## See comment in R/boot.R
    environment(f) <- environment(iofn) <- .GlobalEnv

    if (is.null(cl)) {
      if (is.null(outfile)) {
        outfile <- ""
      }
      cl <- makePSOCKcluster(cores, outfile = outfile)
      on.exit(stopCluster(cl))
    }
    clusterExport(cl,
      varlist = c("fitted", "ytol", "ystep", "trace", "iofn", "n"),
      envir = environment()
    )
    dl <- clusterMap(cl, f, x = xl, h = hl, i = seq_along(xl))
  } else {
    if (!is.null(outfile)) {
      sink(outfile, type = "output")
      sink(outfile, type = "message")
    }
    dl <- switch(parallel,
      multicore = mcmapply(f, x = xl, h = hl, i = seq_along(xl),
                           SIMPLIFY = FALSE, mc.cores = cores),
      serial = Map(f, x = xl, h = hl, i = seq_along(xl))
    )
    if (!is.null(outfile)) {
      sink(type = "output")
      sink(type = "message")
    }
  }

  nrl <- vapply(dl, nrow, 0L)
  out <- data.frame(
    name = rep(factor(profile_names, levels = profile_names), nrl),
    do.call(rbind, dl)
  )
  row.names(out) <- NULL
  class(out) <- c("egf_profile", "data.frame")
  out
}

#' Confidence intervals from likelihood profiles
#'
#' Computes confidence intervals on fixed effect coefficients,
#' log standard deviations of random effects, and linear combinations
#' thereof from their (univariate) likelihood profiles.
#'
#' @param object
#'   An `"egf_profile"` object returned by [profile.egf()].
#' @param parm
#'   Unused argument included for generic consistency.
#' @inheritParams confint.egf
#'
#' @details
#' For each likelihood profile (level of `object$name`),
#' [TMB::confint.tmbprofile()] is called to approximate
#' (via linear interpolation) the two solutions of
#' `deviance(x) = qchisq(level, df = 1)`. These are
#' the lower and upper confidence limits of interest.
#'
#' @return
#' A data frame with character variable `name` equal to
#' `levels(object$name)` and numeric variables `estimate`,
#' `lower`, and `upper` supplying estimates and confidence
#' intervals.
#'
#' @export
#' @importFrom stats qchisq
confint.egf_profile <- function(object, parm, level = 0.95, ...) {
  stop_if_not_in_0_1(level)
  stop_if_not(
    qchisq(level, df = 1) < max(object$deviance, na.rm = TRUE),
    m = paste0(
      "Maximum deviance must exceed `qchisq(level, df = 1)`.\n",
      "Reprofile with higher `max_level` or retry with lower `level`."
    )
  )

  object_split <- split(object, object$name)
  q <- qchisq(level, df = 1)

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
    c(estimate = estimate, lower = lower, upper = upper)
  }
  elu <- do.call(rbind, lapply(object_split, f))

  out <- data.frame(
    name = names(object_split),
    elu,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  attr(out, "level") <- level
  out
}

#' Plot likelihood profiles
#'
#' A method for inspecting computed likelihood profiles.
#'
#' @param x
#'   An `"egf_profile"` object returned by [profile.egf()].
#' @param name
#'   A subset of `levels(x$name)` specifying profiles to be plotted.
#' @param sqrt
#'   A logical scalar. If `TRUE`, then square root-transformed
#'   deviance is plotted.
#' @param level
#'   A numeric vector with elements in (0,1), or otherwise `NULL`.
#'   If non-`NULL` and `sqrt = FALSE`, then line segments are drawn
#'   to show the intersection of the profile with horizontal lines
#'   at `qchisq(level, df = 1)`.
#' @param ...
#'   Optional graphical parameters passed to [graphics::plot()],
#'   such as `type = "o"`. Note that `ann = FALSE` and `axes = FALSE`
#'   are hard-coded. Axes and their labels are drawn separately.
#'
#' @return
#' `NULL` (invisibly).
#'
#' @export
#' @import graphics
#' @importFrom stats confint
plot.egf_profile <- function(x, name = levels(x$name), sqrt = FALSE,
                             level = NULL, ...) {
  stop_if_not(
    is.atomic(name),
    length(name) > 0L,
    name %in% levels(x$name),
    !duplicated(name),
    m = "`name` must be a subset of `levels(x$name)`."
  )
  stop_if_not_tf(sqrt)
  if (!sqrt && !is.null(level)) {
    stop_if_not(
      is.numeric(level),
      length(level) > 0L,
      level > 0,
      level < 1,
      m = paste0(
        "`level` must be NULL or a numeric vector\n",
        "with elements in (0,1)."
      )
    )
  }

  x_split <- split(x, factor(x$name, levels = name))

  f <- if (sqrt) base::sqrt else identity
  ymax <- f(max(x$deviance, na.rm = FALSE))
  ylab <- if (sqrt) expression(sqrt("deviance")) else "deviance"

  if (!sqrt & !is.null(level)) {
    h <- f(qchisq(level, df = 1))
    cil <- lapply(level, function(p) confint(x, level = p))
    m <- match(names(x_split), levels(x$name))
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
      x = x_split[[i]]$value,
      y = f(x_split[[i]]$deviance),
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
    title(xlab = names(x_split)[i], line = 2)
    title(ylab = ylab, line = 2.25)
  }

  invisible(NULL)
}
