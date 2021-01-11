#' Compute likelihood profiles
#'
#' Computes the likelihood profile of specified fixed effect
#' coefficients, log standard deviations of random effects,
#' and linear combinations thereof.
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
    vc <- vcov(fitted)
    hl <- sqrt(diag(vc))[index] / 4
    xl <- unname(index)
  } else {
    ## Covariance matrix of `lin_comb %*% u`
    vc <- lin_comb %*% vcov(fitted) %*% t(lin_comb)
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
    d
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
    name = factor(rep.int(profile_names, nrl), levels = profile_names),
    do.call(rbind, dl)
  )
  row.names(out) <- NULL
  class(out) <- c("egf_profile", "data.frame")
  out
}
