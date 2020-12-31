#' @importFrom TMB MakeADFun
boot <- function(object, n = 6L,
                 method = c("nlminb", "nlm", "Nelder-Mead", "BFGS", "CG"),
                 ...) {
  method <- match.arg(method)
  object$madf_args$parameters <- split(object$par, sub("\\[[0-9]+\\]$", "", names(object$par)))
  if (!has_random(object)) {
    object$madf_args$parameters[c("log_sd", "b")] <- list(NA_real_, NA_real_)
  }
  t(replicate(n, {
    object$madf_args$data$x <- object$madf_out$simulate(object$par)$x
    madf_out_boot <- do.call(MakeADFun, object$madf_args)
    par_name <- switch(method, nlm = "estimate", "par")
    rename_par(optim_madf_out(madf_out_boot, method)[[par_name]])
  }))
}
