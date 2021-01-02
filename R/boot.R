#' @export
#' @importFrom TMB MakeADFun dynlib
#' @importFrom parallel mcmapply makePSOCKcluster clusterExport clusterEvalQ parSapply stopCluster
boot_par <- function(object, n = 6L,
                     method = c("nlminb", "nlm", "Nelder-Mead", "BFGS", "CG"),
                     parallel = c("serial", "multicore", "snow"),
                     cores = getOption("egf.cores", 2L),
                     ...) {
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

  object$tmb_args$parameters <- make_tmb_parameters(
    tmb_data = object$tmb_args$data,
    par_init = object$par
  )

  f <- function(...) {
    tryCatch(
      expr = {
        object$tmb_args$data$x <- object$tmb_out$simulate(object$par)$x
        tmb_out_sim <- do.call(MakeADFun, object$tmb_args)
        optim_tmb_out(tmb_out_sim, method)
        rename_par(tmb_out_sim$env$last.par.best)
      },
      error = function(e) {
        x <- rep(NA_real_, length(object$par))
        names(x) <- names(object$par)
        x
      }
    )
  }
  if (parallel == "multicore") {
    m <- mcmapply(f, seq_len(n), mc.cores = cores)
  } else if (parallel == "snow") {
    dll <- sprintf("%s/epigrowthfit", system.file("src", package = "epigrowthfit"))
    cl <- makePSOCKcluster(cores)
    on.exit(stopCluster(cl))
    var <- c("object", "method", "dll",
             "MakeADFun", "dynlib",
             "optim_tmb_out", "rename_par")
    clusterExport(cl, var, envir = environment())
    clusterEvalQ(cl, dyn.load(dynlib(dll)))
    m <- parSapply(cl, seq_len(n), f)
  } else { # "serial"
    m <- replicate(n, f())
  }
  m
}
