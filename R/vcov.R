#' @export
#' @importFrom TMB sdreport
#' @importFrom stats cov2cor
vcov.egf <- function(object, decontrast = FALSE,
                     cor = FALSE, random = FALSE, ...) {
  stop_if_not_tf(decontrast)
  stop_if_not_tf(cor)
  stop_if_not_tf(random)

  re <- has_random(object)

  if (random) {
    if (!re) {
      stop(
        "`object` does not fit a random effects model.\n",
        "Use `random = FALSE` to retrieve the fixed effects\n",
        "covariance matrix."
      )
    }
    ml <- object$report$cor_list
    rid <- object$tmb_args$data$rid
    names(ml) <- rownames(rid)
    dnl <- lapply(seq_len(nrow(rid)), function(i) {
      pn <- colnames(rid)[rid[i, ] > 0L]
      rep(list(pn), 2L)
    })
    ml <- Map("dimnames<-", ml, dnl)
    if (!cor) {
      sdl <- object$report$sd_list
      names(sdl) <- rownames(rid)
      ml <- Map(function(x, y) outer(x, x) * y, sdl, ml)
    }
    return(ml)
  }

  m <- object$report$cov
  if (decontrast) {
    lc <- make_lin_comb(fnl = object$tmb_args$data$fnl,
                        rid = object$tmb_args$data$rid)
    m <- lc %*% m %*% t(lc)
  }
  if (cor) {
    m <- cov2cor(m)
  }
  dn <- with(object$par_info, {
    c(element[vector == if (decontrast) ".decontrast(beta)" else "beta"],
      if (re) element[vector == "log_sd_b"])
  })
  dimnames(m) <- rep(list(dn), 2L)
  m
}
