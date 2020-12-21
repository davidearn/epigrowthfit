#' @importFrom TMB sdreport
#' @importFrom stats setNames
vcov.egf <- function(object, decontrast = TRUE,
                     full = FALSE, cor = FALSE,
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

    rtl <- rownames(rid)
    if (is.null(random_term_labels)) {
      random_term_labels <- rtl
    } else {
      stop_if_not(
        is.character(random_term_labels),
        length(random_term_labels) > 0L,
        random_term_labels %in% rtl,
        m = paste0(
          "`random_term_labels` must be NULL or a subset\n",
          "of `rownames(object$madf_args$data$rid)`."
        )
      )
    }

    wh <- match(random_term_labels, rtl)
    ml <- object$report$cor_list[wh]
    names(ml) <- rtl[wh]
    dnl <- lapply(wh, function(i) {
      pn <- colnames(rid)[rid[i, ] > 0L]
      rep(list(pn), 2L)
    })
    ml <- Map("dimnames<-", ml, dnl)
    if (!cor) {
      sdl <- object$madf_report$sd_list[wh]
      names(sdl) <- rtl[wh]
      ml <- Map(function(x, y) outer(x, x) * y, sdl, ml)
    }
    return(ml)
  }

  m <- object$report$cov
  dimnames(m) <- rep(list(grep("^b\\[", names(object$par), value = TRUE, invert = TRUE)), 2L)
  k <- grep("^beta\\[", colnames(m))
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
    m <- lc %*% m %*% t(lc) # dimnames discarded
    dimnames(m) <- dn
  }
  if (cor) {
    sd <- sqrt(diag(m))
    m <- m / outer(sd, sd)
  }

  attr(m, "par_info") <- get_par_info(object,
    which = if (full) c("beta", "sd") else "beta",
    decontrast = decontrast
  )
  m
}
