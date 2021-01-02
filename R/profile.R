#' @importFrom tmbprofile
#' @importFrom parallel mcmapply makePSOCKcluster clusterExport clusterEvalQ parSapply stopCluster
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

  if (!is.null(index)) {
    stop_if_not(
      is.numeric(index),
      length(index) > 0L,
      index %in% fitted$nonrandom,
      m = "`index` must be a subset of `fitted$nonrandom`."
    )
    index <- unique(index)
    lin_comb <- NULL
    lin_comb_names <- fitted$par_info$element[index]
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
  } else if (!is.null(parm)) {
    pn <- get_par_names(fitted, link = TRUE)
    stop_if_not(
      is.character(parm),
      length(parm) > 0L,
      (parm <- add_link_string(unique(remove_link_string(parm)))) %in% pn,
      m = "`parm` must be a subset of `get_par_names(fitted)`."
    )
    stop_if_not_tf(decontrast)

    i_parm <- which(object$par_info$response %in% parm)
    i_parm_split <- split(i_parm, object$par_info$vector[i_parm])
    i1 <- unlist(i_parm_split[c("beta", "log_sd_b")])
    i2 <- unlist(i_parm_split[c(".decontrast(beta)", "log_sd_b")])

    if (decontrast) {
      lin_comb <- make_lin_comb(fnl = object$tmb_args$data$fnl,
                                rid = object$tmb_args$data$rid)
      lin_comb <- lin_comb[i1, , drop = FALSE]
      lin_comb_names <- object$par_info$element[i2]
    } else {
      index <- unname(i1)
      lin_comb <- NULL
      lin_comb_names <- object$par_info$element[i1]
    }
  } else {
    stop("One of `index`, `lin_comb`, and `parm` must be non-NULL.")
  }

  if (is.null(lin_comb)) {
    vc <- vcov(fitted)
    h <- sqrt(diag(vc))[index] / 4
    x <- index
  } else {
    vc <- lin_comb %*% vcov(fitted) %*% t(lin_comb)
    h <- sqrt(diag(vc)) / 4
    x <- lapply(seq_len(nrow(lin_comb)), function(i) lin_comb[i, ])
  }
  ytol <- qchisq(max_level, df = 1) / 2
  ystep <- ytol / grid_len

  f <- function(x, h) {
    if (length(x) > 1L) {
      d <- tmbprofile(obj = fitted$tmb_out, lincomb = x, h = h,
                      ytol = ytol, ystep = ystep, trace = trace)
    } else {
      d <- tmbprofile(obj = fitted$tmb_out, name = x, h = h,
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
    on.exit(stopCluster(cl))
    var <- c("fitted", "ytol", "ystep", "trace", "tmbprofile")
    clusterExport(cl, var, envir = environment())
    dl <- clusterMap(cl, f, x = x, h = h)
  } else { # "serial"
    dl <- Map(f, x = x, h = h)
  }

  nr <- vapply(dl, nrow, integer(1L))
  out <- cbind(name = factor(rep(lin_comb_names, nr)), do.call(rbind, dl))
  row.names(out) <- NULL
  class(out) <- c("egf_profile", "data.frame")
  if (!is.null(lin_comb)) {
    attr(out, "lin_comb") <- lin_comb
  }
  out
}
