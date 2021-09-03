#' Parametric bootstrap
#'
#' Sample from the conditional distribution of segments \code{beta} and
#' \code{theta} of the full parameter vector \code{c(beta, theta, b)},
#' given their fitted values.
#'
#' @param object
#'   An \code{"\link{egf}"} object.
#' @param n
#'   A positive integer. The desired number of samples.
#' @param parallel
#'   An \code{"\link{egf_parallel}"} object defining parallelization
#'   options.
#' @param trace
#'   A \link{logical} flag.
#'   If \code{TRUE}, then basic tracing messages indicating progress
#'   are printed. Depending on \code{object$control$trace}, these may
#'   be mixed with optimization output.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A \link{numeric} \link{matrix} inheriting from \link{class}
#' \code{"egf_boot"},
#' with \code{\link{length}(object$best)} rows and \code{n} columns.
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

  tmb_args <- egf_remake_tmb_args(object)
  tmb_out <- object$tmb_out
  best <- object$best
  control <- object$control

  do_sim <- function(i) {
    if (trace) {
      cat(sprintf("Commencing bootstrap simulation %d of %d...\n", i, n))
    }
    on <- openmp(n = NULL)
    if (on > 0L) {
      openmp(n = control$omp_num_threads)
      on.exit(openmp(n = on))
    }
    ## Update
    tmb_args$data$x <- tmb_out$simulate(best)$x
    ## Retape
    tmb_out_retape <- do.call(MakeADFun, tmb_args)
    ## Patch
    tmb_out_retape$fn <- egf_patch_fn(tmb_out_retape$fn, inner_optimizer = control$inner_optimizer)
    tmb_out_retape$gr <- egf_patch_gr(tmb_out_retape$gr, inner_optimizer = control$inner_optimizer)
    ## Optimize
    optimizer <- control$optimizer$f
    optimizer_args <- c(
      tmb_out[c("par", "fn", "gr")],
      control$optimizer["control"],
      control$optimizer[["args"]]
    )
    tryCatch(
      expr = {
        do.call(optimizer, optimizer_args)
        tmb_out_retape$env$last.par.best
      },
      error = function(cond) {
        cat(sprintf("Error in bootstrap simulation %d of %d:\n %s", i, n, conditionMessage(cond)))
        rep_len(NaN, length(best))
      }
    )
  }

  if (parallel$method == "snow") {
    ## Doing this to avoid duplication of 'clusterExport'
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
    clusterExport(cl,
      varlist = c(
        "tmb_args", "trace", "n",
        ## Could simply serialize 'object' but trying to minimize overhead
        "tmb_out", "best", "control",
        ## 'egf_patch_fn' and 'egf_patch_gr' are internal functions, so they
        ## must be made available in global environment of worker processes
        "egf_patch_fn", "egf_patch_gr"
      ),
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
      serial = vapply(seq_len(n), do_sim, unname(best))
    )
    if (nzchar(parallel$outfile)) {
      sink(type = "output")
      sink(type = "message")
    }
  }
  rownames(res) <- names(best)
  class(res) <- c("egf_boot", "matrix", "array")
  res
}
