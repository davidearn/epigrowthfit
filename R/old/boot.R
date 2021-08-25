#' Parametric bootstrap
#'
#' Sample from the conditional distribution of the full parameter vector
#' given its fitted value.
#'
#' @param object
#'   An \code{"\link{egf}"} object.
#' @param n
#'   A positive integer. The desired number of samples.
#' @param parallel
#'   An \code{"\link{egf_parallel}"} object defining parallelization options.
#' @param trace
#'   A \link{logical} flag. If \code{TRUE}, then basic tracing messages
#'   indicating progress are printed. Depending on \code{object$control$trace},
#'   these may be mixed with optimization output.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A \link{numeric} \link{matrix} inheriting from \link{class}
#' \code{"egf_boot"}, with \code{\link{length}(object$best)} rows
#' and \code{n} columns.
#'
#' @export
#' @importFrom TMB MakeADFun
#' @import parallel
#' @importFrom TMB openmp
boot_par <- function(object,
                     n = 6L,
                     parallel = egf_parallel(),
                     trace = FALSE,
                     ...) {
  stopifnot(
    inherits(object, "egf"),
    inherits(parallel, "egf_parallel")
  )
  stop_if_not_integer(n, "positive")
  stop_if_not_true_false(trace)

  ## Make sure that bootstrap optimizations start from fitted model
  object$tmb_args$parameters <- egf_make_tmb_parameters(
    tmb_data = object$tmb_args$data,
    fit = TRUE,
    init = object$best
  )

  do_sim <- function(i) {
    if (trace) {
      cat(sprintf("Commencing bootstrap simulation %d of %d...\n", i, n))
    }
    on <- openmp(n = NULL)
    if (on > 0L) {
      openmp(n = object$control$omp_num_threads)
      on.exit(openmp(n = on))
    }
    ## Simulate data from fitted model
    object$tmb_args$data$x <- object$tmb_out$simulate(object$best)$x
    ## Define objective function and gradient conditional on simulated data
    tmb_out_sim <- do.call(MakeADFun, object$tmb_args)
    ## Patch
    tmb_out_sim$fn <- egf_patch_fn(tmb_out_sim$fn, inner_optimizer = object$control$inner_optimizer)
    tmb_out_sim$gr <- egf_patch_gr(tmb_out_sim$gr, inner_optimizer = object$control$inner_optimizer)
    ## Optimize
    optimizer <- object$control$optimizer$f
    optimizer_args <- c(
      tmb_out_sim[c("par", "fn", "gr")],
      object$control$optimizer["control"],
      object$control$optimizer[["args"]]
    )
    tryCatch(
      expr = {
        do.call(optimizer, optimizer_args)
        tmb_out_sim$env$last.par.best
      },
      error = function(cond) {
        cat(sprintf("Error in bootstrap simulation %d of %d:\n %s", i, n, conditionMessage(cond)))
        rep_len(NaN, length(object$best))
      }
    )
  }

  if (parallel$method == "snow") {
    ## Doing this to avoid duplication of `clusterExport`
    ## (i.e., unnecessary serialization of function environment):
    ## https://stackoverflow.com/questions/17402077/how-to-clusterexport-a-function-without-its-evaluation-environment
    ## https://stackoverflow.com/questions/18035711/environment-and-scope-when-using-parallel-functions
    environment(do_sim) <- .GlobalEnv

    if (is.null(parallel$cl)) {
      cl <- do.call(makePSOCKcluster, parallel$args)
      on.exit(stopCluster(cl))
    } else {
      cl <- parallel$cl
    }
    clusterEvalQ(cl, library("TMB"))
    ## Need to make internal functions `egf_patch_fn` and `egf_patch_gr`
    ## available in global environment of worker processes
    clusterExport(cl,
      varlist = c("object", "n", "trace", "egf_patch_fn", "egf_patch_gr"),
      envir = environment()
    )
    clusterSetRNGStream(cl)
    res <- parSapply(cl, seq_len(n), do_sim)
  } else {
    if (nzchar(parallel$outfile)) {
      outfile <- file(parallel$outfile, open = "wt")
      sink(outfile, type = "output")
      sink(outfile, type = "message")
    }
    res <- switch(parallel$method,
      multicore = simplify2array(do.call(mclapply, c(list(X = seq_len(n), FUN = do_sim), parallel$args))),
      serial = vapply(seq_len(n), do_sim, unname(object$best))
    )
    if (nzchar(parallel$outfile)) {
      sink(type = "output")
      sink(type = "message")
    }
  }
  rownames(res) <- names(object$best)
  class(res) <- c("egf_boot", "matrix", "array")
  res
}


