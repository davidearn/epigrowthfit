#' Parametric bootstrap
#'
#' Sample from the conditional distribution of the full parameter vector
#' given the fitted epidemic model.
#'
#' @inheritParams egf
#' @param object
#'   An `"egf"` object returned by [egf()].
#' @param n
#'   A positive integer. The desired number of samples.
#' @param trace
#'   A logical scalar. If `TRUE`, then basic tracing information is
#'   printed.
#' @param parallel
#'   A character string indicating a library for parallel computation.
#'   `"serial"` indicates no parallelization. `"multicore"` forks and
#'   is intended for use from a terminal rather than, e.g., RStudio.
#'   On Windows, it is equivalent to `"serial"`. `"snow"` creates
#'   socket clusters and parallelizes on both Unix-alikes and Windows.
#' @param cores
#'   A positive integer. The number of worker processes spawned when
#'   `parallel != "serial"`. See also [parallel::detectCores()].
#' @param outfile
#'   A character string containing a file path or otherwise `NULL`.
#'   If a file path, then console output is diverted there. If `NULL`,
#'   then there is no diversion. When `parallel = "snow"`, diversion
#'   may be necessary to view output.
#' @param cl
#'   A optional socket cluster created by
#'   [parallel::makePSOCKcluster()], to be used
#'   when `parallel = "snow"`. If non-`NULL`,
#'   then `cores` and `outfile` are ignored.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A numeric matrix inheriting from class `"egf_boot"`,
#' with `length(object$par)` rows and `n` columns.
#'
#' @export
#' @import parallel
boot_par <- function(object, n = 6L,
                     method = c("nlminb", "nlm", "Nelder-Mead", "BFGS"),
                     trace = TRUE,
                     parallel = c("serial", "multicore", "snow"),
                     cores = getOption("egf.cores", 2L),
                     outfile = NULL,
                     cl = NULL,
                     ...) {
  ## FIXME: Should this be a method for a boot() generic?
  stop_if_not(
    inherits(object, "egf"),
    m = "`object` must be an \"egf\" object."
  )
  stop_if_not_positive_integer(n)
  method <- match.arg(method)
  stop_if_not_tf(trace)
  parallel <- match.arg(parallel)
  check_parallel(parallel, cores, outfile, cl)

  ## Make sure that bootstrap optimizations start from
  ## the fitted parameter vector
  object$tmb_args$parameters <- make_tmb_parameters(
    tmb_data = object$tmb_args$data,
    par_init = object$par
  )

  iofn <- function(i) sprintf("%*d of %d", nchar(n), i, n)
  f <- function(i) {
    if (trace) cat("Starting bootstrap simulation", iofn(i), "...\n")
    ## Simulate new data from full parameter vector
    object$tmb_args$data$x <- object$tmb_out$simulate(object$par)$x
    ## Redefine objective function
    tmb_out_sim <- do.call(MakeADFun, object$tmb_args)
    ## Optimize without stopping in the event of NaN
    ## function evaluations or similar errors
    tryCatch(
      expr = {
        optim_tmb_out(tmb_out_sim, method)
        tmb_out_sim$env$last.par.best
      },
      error = function(e) {
        cat("Error in bootstrap simulation", iofn(i), ":\n", e[[1L]])
        rep.int(NA_real_, length(object$par))
      }
    )
  }

  if (parallel == "snow") {
    ## Doing this to avoid duplication of clusterExport().
    ## Function environments would otherwise be serialized
    ## unnecessarily.
    environment(f) <- environment(iofn) <- .GlobalEnv

    if (is.null(cl)) {
      if (is.null(outfile)) {
        outfile <- ""
      }
      cl <- makePSOCKcluster(cores, outfile = outfile)
      on.exit(stopCluster(cl))
    }
    clusterEvalQ(cl, library("epigrowthfit"))
    ## Need to make internal function optim_tmb_out()
    ## available in global environment of workers
    clusterExport(cl,
      varlist = c("object", "method", "trace", "iofn", "n", "optim_tmb_out"),
      envir = environment()
    )
    clusterSetRNGStream(cl)
    m <- parSapply(cl, seq_len(n), f)
  } else {
    if (!is.null(outfile)) {
      sink(outfile, type = "output")
      sink(outfile, type = "message")
    }
    m <- switch(parallel,
      multicore = simplify2array(mclapply(seq_len(n), f, mc.cores = cores)),
      serial = vapply(seq_len(n), f, numeric(length(object$par)))
    )
    if (!is.null(outfile)) {
      sink(type = "output")
      sink(type = "message")
    }
  }
  rownames(m) <- names(object$par)
  class(m) <- c("egf_boot", "matrix", "array")
  m
}
