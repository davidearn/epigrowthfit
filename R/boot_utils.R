#' Check parallelization options
#'
#' Validate arguments used by most methods supporting parallelization.
#'
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
#'
#' @return
#' `NULL` (invisibly).
#'
#' @keywords internal
check_parallel <- function(parallel, cores, outfile, cl) {
  if (parallel == "multicore" || (parallel == "snow" && is.null(cl))) {
    stop_if_not_positive_integer(cores, n = 2L)
    if (!is.null(outfile)) {
      stop_if_not_character_string(outfile, n = 2L)
    }
  } else if (parallel == "snow" && !is.null(cl)) {
    stop_if_not(
      inherits(cl, "SOCKcluster"),
      m = "`cl` must be a socket cluster or NULL.",
      n = 2L
    )
  }
  invisible(NULL)
}
