#' @export
#' @importFrom parallel mcmapply makePSOCKcluster clusterExport clusterEvalQ parSapply stopCluster
boot_par <- function(object, n = 6L,
                     method = c("nlminb", "nlm", "Nelder-Mead", "BFGS", "CG"),
                     parallel = c("serial", "multicore", "snow"),
                     cores = getOption("egf.cores", 2L),
                     ...) {
  ## FIXME: Export a boot() generic from another package?
  stop_if_not(
    inherits(object, "egf"),
    m = "`object` must be an \"egf\" object."
  )
  stop_if_not_positive_integer(n)
  method <- match.arg(method)
  parallel <- match.arg(parallel)
  if (parallel != "serial") {
    stop_if_not_positive_integer(cores)
  }

  ## Make sure that bootstrap optimizations start from
  ## fitted parameter vector
  object$tmb_args$parameters <- make_tmb_parameters(
    tmb_data = object$tmb_args$data,
    par_init = object$par
  )

  f <- function(...) {
    ## Simulate new data from fitted parameter vector
    object$tmb_args$data$x <- object$tmb_out$simulate(object$par)$x
    ## Redefine objective function
    tmb_out_sim <- do.call(TMB::MakeADFun, object$tmb_args)
    ## Optimize without stopping in the event of NaN
    ## function evaluations or similar errors
    tryCatch(
      expr = {
        epigrowthfit:::optim_tmb_out(tmb_out_sim, method)
        epigrowthfit:::rename_par(tmb_out_sim$env$last.par.best)
      },
      error = function(e) {
        x <- rep(NA_real_, length(object$par))
        names(x) <- names(object$par)
        x
      }
    )
  }

  if (parallel == "multicore") {
    mcmapply(f, seq_len(n), mc.cores = cores)
  } else if (parallel == "snow") {
    dll <- sprintf("%s/epigrowthfit", system.file("src", package = "epigrowthfit"))
    cl <- makePSOCKcluster(cores)
    on.exit(stopCluster(cl))
    var <- c("object", "method", "dll")
    clusterExport(cl, var, envir = environment())
    clusterEvalQ(cl, dyn.load(TMB::dynlib(dll)))
    parSapply(cl, seq_len(n), f)
  } else { # "serial"
    replicate(n, f())
  }
}
