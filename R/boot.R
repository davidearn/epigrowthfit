#' Parametric bootstrap
#'
#' Sample from the conditional distribution of segments \code{beta} and
#' \code{theta} of the full parameter vector \code{c(beta, theta, b)}
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
#' A \link[=double]{numeric} \link{matrix} inheriting from \link{class}
#' \code{"egf_boot"}, with \code{\link{length}(object$best)} rows and
#' \code{n} columns.
#'
#' @examples
#' example("egf", "epigrowthfit")
#' set.seed(181952L)
#' bb <- boot_par(object, n = 6L)
#' matplot(t(bb[object$nonrandom, ]), type = "o", las = 1,
#'   xlab = "simulation",
#'   ylab = "parameter value"
#' )
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

  set_omp_num_threads <- quote({
    current <- openmp(n = NULL)
    if (current > 0L) {
      openmp(n = object$control$omp_num_threads)
      on.exit(openmp(n = current), add = TRUE)
    }
  })
  eval(set_omp_num_threads)

  ## Simulations are fast and can be done in the master process to avoid
  ## serializing a TMB object
  X <- replicate(n, object$tmb_out$simulate(object$best)$x, simplify = FALSE)

  ## Reconstruct list of arguments to 'MakeADFun' from object internals
  ## for retaping
  tmb_args <- egf_remake_tmb_args(object)

  do_sim <- function(i, x) {
    if (trace) {
      cat(sprintf("Commencing bootstrap simulation %d of %d...\n", i, n))
    }
    ## Update
    tmb_args$data$x <- x
    ## Retape
    tmb_out_retape <- do.call(MakeADFun, tmb_args)
    ## Patch
    tmb_out_retape$fn <- egf_patch_fn(tmb_out_retape$fn, inner_optimizer = object$control$inner_optimizer)
    tmb_out_retape$gr <- egf_patch_gr(tmb_out_retape$gr, inner_optimizer = object$control$inner_optimizer)
    ## Optimize
    optimizer <- object$control$optimizer$f
    optimizer_args <- c(
      tmb_out_retape[c("par", "fn", "gr")],
      object$control$optimizer["control"],
      object$control$optimizer[["args"]]
    )
    tryCatch(
      expr = {
        do.call(optimizer, optimizer_args)
        tmb_out_retape$env$last.par.best
      },
      error = function(cond) {
        cat(sprintf("Error in bootstrap simulation %d of %d:\n %s", i, n, conditionMessage(cond)))
        rep_len(NaN, length(object$best))
      }
    )
  }

  if (parallel$method == "snow") {
    ## We use 'clusterExport' to export necessary objects to the
    ## global environments of all worker processes. As a result,
    ## function environments are unused and need not be serialized.
    ## By replacing them with the global environment,
    ## which is never serialized, we avoid unnecessary overhead.
    ## https://stackoverflow.com/questions/17402077/how-to-clusterexport-a-function-without-its-evaluation-environment
    ## https://stackoverflow.com/questions/18035711/environment-and-scope-when-using-parallel-functions
    environment(do_sim) <- .GlobalEnv

    ## Only these elements of 'object' are necessary
    object <- object[c("control", "best")]

    if (is.null(parallel$cl)) {
      cl <- do.call(makePSOCKcluster, parallel$args)
      on.exit(stopCluster(cl), add = TRUE)
    } else {
      cl <- parallel$cl
    }
    clusterExport(cl,
      varlist = c("tmb_args", "object", "trace", "n", "set_omp_num_threads"),
      envir = environment()
    )
    clusterEvalQ(cl, {
      library("TMB")
      eval(set_omp_num_threads)
      ## Hack to avoid 'check' note about 'epigrowthfit:::<internal function>'
      egf_patch_fn <- eval(call(":::", quote(epigrowthfit), quote(egf_patch_fn)))
      egf_patch_gr <- eval(call(":::", quote(epigrowthfit), quote(egf_patch_gr)))
    })
    clusterSetRNGStream(cl)
    res <- clusterMap(cl, do_sim, i = seq_len(n), x = X, simplify = TRUE)
  } else {
    if (nzchar(parallel$outfile)) {
      outfile <- file(parallel$outfile, open = "wt")
      sink(outfile, type = "output")
      sink(outfile, type = "message")
      on.exit(add = TRUE, {
        sink(type = "output")
        sink(type = "message")
      })
    }
    res <- switch(parallel$method,
      multicore = do.call(mcmapply, c(list(FUN = do_sim, i = seq_len(n), x = X), parallel$args)),
      serial = mapply(do_sim, i = seq_len(n), x = X)
    )
  }
  rownames(res) <- names(object$best)
  class(res) <- c("egf_boot", "matrix", "array")
  res
}
