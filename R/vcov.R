#' @export
#' @importFrom TMB sdreport
#' @importFrom stats cov2cor
vcov.egf <- function(object, decontrast = FALSE,
                     cor = FALSE, random = FALSE, ...) {
  stop_if_not_tf(decontrast)
  stop_if_not_tf(cor)
  stop_if_not_tf(random)

  if (random) {
    if (!has_random(object)) {
      stop(
        "`object` does not fit a random effects model.\n",
        "Use `random = FALSE` to retrieve the fixed effects\n",
        "covariance matrix."
      )
    }

    ## List of random effects correlation matrices. There is one
    ## matrix for each unique random effects term in the `p` mixed
    ## effects models. The number of columns in a given matrix
    ## is equal to the number of model formulae (out of `p`) that
    ## specify the corresponding term as a random effect.
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

  ## Covariance matrix of `u = c(beta, log_sd_b)`
  m <- object$report$cov

  if (decontrast) {
    ## Matrix of linear coefficients such that `lc %*% u`
    ## returns `v = c(decontrast(beta), log_sd_b)`
    lc <- make_lin_comb(fnl = object$tmb_args$data$fnl,
                        rid = object$tmb_args$data$rid)
    ## Covariance matrix of `v`
    m <- lc %*% m %*% t(lc)
  }

  if (cor) {
    m <- cov2cor(m)
  }

  dn <- with(object$par_info, {
    c(name[vector == if (decontrast) ".decontrast(beta)" else "beta"],
      name[vector == "log_sd_b"])
  })
  dimnames(m) <- rep(list(dn), 2L)
  m
}
