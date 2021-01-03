#' @export
#' @importFrom parallel mcmapply makePSOCKcluster clusterExport clusterMap stopCluster
profile.egf <- function(fitted, index = NULL, lin_comb = NULL,
                        parm = get_par_names(fitted), decontrast = FALSE,
                        max_level = 0.99, grid_len = 12,
                        trace = FALSE,
                        parallel = c("serial", "multicore", "snow"),
                        cores = getOption("egf.cores", 2L),
                        ...) {
  stop_if_not_tf(trace)
  stop_if_not(
    is.numeric(max_level),
    length(max_level) == 1L,
    is.finite(max_level),
    max_level > 0 && max_level < 1,
    m = "`max_level` must be a number in the interval (0,1)."
  )
  stop_if_not(
    is.numeric(grid_len),
    length(grid_len) == 1L,
    is.finite(grid_len),
    grid_len >= 1,
    m = "`grid_len` must be 1 or greater."
  )
  parallel <- match.arg(parallel)
  if (parallel != "serial") {
    stop_if_not_positive_integer(cores)
  }

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
    lin_comb_names <- fitted$par_info$name[index]

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
    lin_comb_names <- seq_len(nrow(lin_comb))

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

    ## If profiling elements of `v = c(decontrast(beta), log_sd_b)`
    if (decontrast) {
      ## Matrix of linear coefficients such that
      ## `lin_comb %*% u = v[relevant]`
      lin_comb <- make_lin_comb_for_parm(object, parm)
      lin_comb_names <- rownames(lin_comb)
    ## If profiling elements of `u`
    } else {
      index <- make_index_for_parm(object, parm)
      lin_comb <- NULL
      lin_comb_names <- names(index)
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

  f <- function(x, h) {
    if (length(x) > 1L) {
      d <- TMB::tmbprofile(obj = fitted$tmb_out, lincomb = x, h = h,
                           ytol = ytol, ystep = ystep, trace = trace)
    } else {
      d <- TMB::tmbprofile(obj = fitted$tmb_out, name = x, h = h,
                           ytol = ytol, ystep = ystep, trace = trace)
    }
    d[[2L]] <- 2 * (d[[2L]] - min(d[[2L]], na.rm = TRUE)) # deviance = 2 * diff(nll)
    names(d) <- c("value", "deviance")
    d
  }

  if (parallel == "multicore") {
    dl <- mcmapply(f, x = xl, h = hl, SIMPLIFY = FALSE, mc.cores = cores)
  } else if (parallel == "snow") {
    cl <- makePSOCKcluster(cores)
    on.exit(stopCluster(cl))
    var <- c("fitted", "ytol", "ystep", "trace")
    clusterExport(cl, var, envir = environment())
    dl <- clusterMap(cl, f, x = xl, h = hl)
  } else { # "serial"
    dl <- Map(f, x = xl, h = hl)
  }

  nrl <- vapply(dl, nrow, integer(1L))
  out <- data.frame(name = factor(rep(lin_comb_names, nrl)), do.call(rbind, dl))
  row.names(out) <- NULL
  class(out) <- c("egf_profile", "data.frame")
  out
}
