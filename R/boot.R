#' Parametric bootstrap
#'
#' Sample from the conditional distribution of the full
#' parameter vector given the fitted epidemic model.
#'
#' @inheritParams egf
#' @param object
#'   An `"egf"` object returned by [egf()].
#' @param n
#'   A positive integer. The desired number of samples.
#' @param trace
#'   A logical scalar. If `TRUE`, then basic tracing messages
#'   are printed.
#' @inheritParams check_parallel
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A numeric matrix inheriting from class `"egf_boot"`,
#' with `length(object$best)` rows and `n` columns.
#'
#' @export
#' @import parallel
boot_par <- function(object,
                     n = 6L,
                     method = c("nlminb", "nlm", "Nelder-Mead", "BFGS"),
                     trace = TRUE,
                     parallel = c("serial", "multicore", "snow"),
                     cores = getOption("egf.cores", 2L),
                     outfile = NULL,
                     cl = NULL,
                     ...) {
  ## FIXME: Should this be a method for a `boot()` generic?
  stop_if_not(
    inherits(object, "egf"),
    m = "`object` must be an \"egf\" object."
  )
  stop_if_not_positive_integer(n)
  method <- match.arg(method)
  stop_if_not_true_false(trace)
  parallel <- match.arg(parallel)
  check_parallel(parallel, cores, outfile, cl)

  ## Make sure that bootstrap optimizations start from fitted model
  object$tmb_args$parameters <- make_tmb_parameters(
    tmb_data = object$tmb_args$data,
    init = object$best
  )

  i_of_n <- function(i) {
    sprintf("%*d of %d", nchar(n), i, n)
  }
  do_sim <- function(i) {
    if (trace) {
      cat("Starting bootstrap simulation", i_of_n(i), "...\n")
    }
    ## Simulate new data from fitted model
    object$tmb_args$data$x <- object$tmb_out$simulate(object$best)$x
    ## Define objective function given simulated data
    tmb_out_sim <- do.call(MakeADFun, object$tmb_args)
    ## Don't stop in the event of NaN function evaluations
    ## or similar errors
    tryCatch(
      expr = {
        optim_tmb_out(tmb_out_sim, method)
        tmb_out_sim$env$last.par.best
      },
      error = function(e) {
        cat("Error in bootstrap simulation", i_of_n(i), ":\n", e[[1L]])
        rep_len(NA_real_, length(object$best))
      }
    )
  }

  if (parallel == "snow") {
    ## Doing this to avoid duplication of `clusterExport()`.
    ## Function environments would otherwise be serialized
    ## unnecessarily.
    environment(do_sim) <- environment(i_of_n) <- .GlobalEnv

    if (is.null(cl)) {
      if (is.null(outfile)) {
        outfile <- ""
      }
      cl <- makePSOCKcluster(cores, outfile = outfile)
      on.exit(stopCluster(cl))
    }
    clusterEvalQ(cl, library("epigrowthfit"))
    ## Need to make internal function `optim_tmb_out()`
    ## available in global environment of worker processes
    clusterExport(cl,
      varlist = c("object", "method", "trace", "i_of_n", "n", "optim_tmb_out"),
      envir = environment()
    )
    clusterSetRNGStream(cl)
    m <- parSapply(cl, seq_len(n), do_sim)
  } else {
    if (!is.null(outfile)) {
      sink(outfile, type = "output")
      sink(outfile, type = "message")
    }
    m <- switch(parallel,
      multicore = simplify2array(mclapply(seq_len(n), do_sim, mc.cores = cores)),
      serial = vapply(seq_len(n), do_sim, rep_len(0, length(object$best)))
    )
    if (!is.null(outfile)) {
      sink(type = "output")
      sink(type = "message")
    }
  }
  rownames(m) <- names(object$best)
  class(m) <- c("egf_boot", "matrix", "array")
  m
}
