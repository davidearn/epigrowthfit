#' @importFrom TMB sdreport
#' @importFrom stats setNames
vcov.egf <- function(object, decontrast = TRUE, full = FALSE, cor = FALSE,
                     random = FALSE, random_term_labels = NULL, ...) {
  for (s in c("decontrast", "full", "cor", "random")) {
    a <- get(s, inherits = FALSE)
    stop_if_not(
      is.logical(a),
      length(a) == 1L,
      !is.na(a),
      m = sprintf("`%s` must be TRUE or FALSE.", s)
    )
  }

  re <- has_random(object)
  rid <- object$madf_args$data$rid

  if (random) {
    if (!re) {
      stop(
        "`object` does not fit a random effects model.\n",
        "Use `random = FALSE` to retrieve the fixed effects\n",
        "variance-covariance matrix."
      )
    }

    rtl0 <- rownames(rid)
    if (is.null(random_term_labels)) {
      random_term_labels <- rtl0
    } else {
      stop_if_not(
        is.character(random_term_labels),
        length(random_term_labels) > 0L,
        random_term_labels %in% rtl0,
        m = paste0(
          "`random_term_labels` must be a subset of\n",
          "`rownames(object$madf_args$data$rid)`."
        )
      )
    }

    wh <- match(random_term_labels, rtl0)
    ml <- object$madf_report$cor_list[wh]
    names(ml) <- rtl0[wh]
    dnl <- lapply(wh, function(i) {
      pn <- colnames(rid)[rid[i, ] > 0L]
      rep(list(pn), 2L)
    })
    ml <- Map("dimnames<-", ml, dnl)
    if (cor) {
      return(ml)
    }
    sdl <- object$madf_report$sd_list[wh]
    names(sdl) <- rtl0[wh]
    return(Map(function(x, y) outer(x, x) * y, sdl, ml))
  }

  m <- object$madf_sdreport$cov.fixed
  k <- which(colnames(m) == "beta")
  if (re && !full) {
    m <- m[k, k]
  }
  if (decontrast) {
    lc <- matrix(0L, nrow = nrow(m), ncol = ncol(m))
    j <- 1L
    for (nl in object$madf_args$data$fnl) {
      i <- seq.int(from = j, length.out = nl)
      lc[i, j] <- 1L
      if (nl > 1L) {
        lc[i, i[-1L]] <- contr.sum(nl)
      }
      j <- j + nl
    }
    if (re && full) {
      lc[-k, -k] <- diag(rep(1L, nrow(m) - length(k)))
    }
    dn <- dimnames(m)
    m <- lc %*% m %*% t(lc)
    dimnames(m) <- dn
  }
  if (cor) {
    sd <- sqrt(diag(m))
    m <- m / outer(sd, sd)
  }

  attr(m, "par_info") <- get_par_info(
    par = object$par,
    madf_data = object$madf_args$data,
    which = if (full) c("beta", "sd_b") else "beta",
    decontrast = decontrast
  )
  m
}
