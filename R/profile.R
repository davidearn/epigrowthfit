#' @importFrom parallel mcmapply makePSOCKcluster clusterMap stopCluster
profile.egf <- function(fitted, index = NULL, lin_comb = NULL,
                        parm = "r", decontrast = TRUE,
                        max_level = 0.99, grid_len = 12,
                        trace = FALSE,
                        parallel = c("serial", "multicore", "snow"),
                        cores = getOption("egf.cores", 2L),
                        ...) {
  stop_if_not(
    is.numeric(trace) || is.logical(trace),
    length(trace) == 1L,
    !is.na(trace),
    m = "`trace` must be a numeric or logical scalar."
  )
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
    stop_if_not(
      is.numeric(cores),
      length(cores) == 1L,
      is.finite(cores),
      cores >= 1,
      m = "`cores` must be 1 or greater."
    )
  }

  n <- length(fitted$inr)

  if (!is.null(index)) {
    stop_if_not(
      is.numeric(index),
      length(index) > 0L,
      index %in% fitted$inr,
      m = "`index` must be a subset of `fitted$inr`."
    )
    index <- unique(index)
    pin <- droplevels(fitted$par_info[index, -1L])
    lin_comb <- NULL
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
      ncol(lin_comb) == n,
      apply(lin_comb, 1L, function(x) any(x != 0)),
      m = sprintf("`lin_comb` must have %d columns and at least\none nonzero element per row.", n)
    )
    pin <- NULL
  } else if (!is.null(parm)) {
    stop_if_not(
      is.logical(decontrast),
      length(decontrast) == 1L,
      !is.na(decontrast),
      m = "`decontrast` must be TRUE or FALSE."
    )
    pn <- sub("^log_", "", get_par_names(fitted))
    stop_if_not(
      is.character(parm),
      length(parm) > 0L,
      parm %in% pn,
      m = paste0(
        "`parm` must be a subset of:\n",
        paste(sprintf("\"%s\"", pn), collapse = ", ")
      )
    )
    parm <- unique(parm)
    parm <- sprintf("log_%s", parm)

    ## Matrix of linear coefficients if `decontrast = TRUE`,
    ## index vector if `decontrast = FALSE`
    lin_comb <- make_lin_comb(fitted, parm, decontrast)
    if (is.vector(lin_comb)) {
      index <- lin_comb
      lin_comb <- NULL
      pin <- droplevels(fitted$par_info[index, -1L])
    } else {
      a <- attributes(lin_comb)
      pin <- droplevels(a$par_info[a$index, -1L])
    }
  } else {
    stop("One of `index`, `lin_comb`, and `parm` must be non-NULL.")
  }

  if (is.null(lin_comb)) {
    vc <- vcov(fitted, full = TRUE)
    h <- sqrt(diag(vc))[index] / 4
    x <- index
  } else {
    vc <- lin_comb %*% vcov(fitted, full = TRUE) %*% t(lin_comb)
    h <- sqrt(diag(vc)) / 4
    x <- lapply(seq_len(nrow(lin_comb)), function(i) lin_comb[i, ])
  }
  ytol <- qchisq(max_level, df = 1) / 2
  ystep <- ytol / grid_len

  f <- function(x, h) {
    if (length(x) > 1L) {
      d <- TMB::tmbprofile(obj = fitted$madf_out, lincomb = x, h = h,
                           ytol = ytol, ystep = ystep, trace = trace)
    } else {
      d <- TMB::tmbprofile(obj = fitted$madf_out, name = x, h = h,
                           ytol = ytol, ystep = ystep, trace = trace)
    }
    d[[2L]] <- 2 * (d[[2L]] - min(d[[2L]], na.rm = TRUE))
    names(d) <- c("value", "deviance")
    d
  }

  if (parallel == "multicore") {
    dl <- mcmapply(f, x = x, h = h, SIMPLIFY = FALSE, mc.cores = cores)
  } else if (parallel == "snow") {
    cl <- makePSOCKcluster(cores)
    dl <- clusterMap(cl, f, x = x, h = h)
    stopCluster(cl)
  } else { # "serial"
    dl <- Map(f, x = x, h = h)
  }

  nr <- vapply(dl, nrow, integer(1L))
  ## If user gave `index`    or `parm` and `decontrast = TRUE`
  if (is.null(lin_comb)) {
    index_long <- rep(index, nr)
  ## Otherwise
  ## If user gave `lin_comb` or `parm` and `decontrast = FALSE`
  } else {
    index_long <- rep(seq_len(nrow(lin_comb)), nr) # dummy index
  }

  ## If user gave `lin_comb`
  if (is.null(pin)) {
    out <- cbind(index = index_long, do.call(rbind, dl))
  ## Otherwise
  } else {
    pin_long <- pin[rep(seq_len(nrow(pin)), nr), ]
    out <- cbind(index = index_long, pin_long, do.call(rbind, dl))
  }

  row.names(out) <- NULL
  class(out) <- c("egf_profile", "data.frame")
  if (!is.null(lin_comb)) {
    attr(out, "lin_comb") <- lin_comb
  }
  out
}
