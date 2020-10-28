#' Print models of epidemic growth
#'
#' @description
#' Methods for printing objects of class "egf_init" or "egf".
#'
#' @param x An "egf_init" or "egf" object.
#' @param ... Unused optional arguments.
#'
#' @return
#' `x` (invisibly).
#'
#' @seealso [egf_init()], [egf()]
#' @name print.egf
NULL

#' @rdname print.egf
#' @export
print.egf_init <- function(x, ...) {
  cstr <- switch(x$curve,
    exponential = "an exponential model",
    logistic    = "a logistic model",
    richards    = "a Richards model"
  )
  bstr <- if (x$include_baseline) {
    "with a linear baseline and"
  } else {
    "with"
  }
  dstr <- switch(x$distr,
    poisson = "Poisson-distributed observations.",
    nbinom  = "negative binomial observations."
  )
  uvec <- c(r = "per day", thalf = "days", b = "per day")
  uvec <- uvec[names(uvec) %in% names(x$theta_init)]
  cat("Pass this \"egf_init\" object to `egf()` to fit", cstr, "\n")
  cat(bstr, dstr, "\n")
  cat("\n")
  cat("Fitting window:\n")
  cat("\n")
  cat("index   ", x$first, ":", x$last, "\n", sep = "")
  cat(" date   (", as.character(x$date[x$first]), ", ", as.character(x$date[x$last+1]), "]\n", sep = "")
  cat("cases   ", sum(x$cases[x$first:x$last]), " of ", sum(x$cases), "\n", sep = "")
  cat("\n")
  cat("Initial parameter estimates:\n")
  cat("\n")
  print(x$theta_init)
  cat("\n")
  cat("Units:\n")
  cat("\n")
  print(uvec, quote = FALSE)
  invisible(x)
}

#' @rdname print.egf
#' @export
print.egf <- function(x, ...) {
  y <- x$init
  cstr <- switch(y$curve,
    exponential = "an exponential model",
    logistic    = "a logistic model",
    richards    = "a Richards model"
  )
  bstr <- if (y$include_baseline) {
    "with a linear baseline and"
  } else {
    "with"
  }
  dstr <- switch(y$distr,
    poisson = "Poisson-distributed observations.",
    nbinom  = "negative binomial observations."
  )
  uvec <- c(r = "per day", thalf = "days", b = "per day")
  uvec <- uvec[names(uvec) %in% names(x$theta_fit)]
  cat("This \"egf\" object fits", cstr, "\n")
  cat(bstr, dstr, "\n")
  cat("\n")
  cat("Fitting window:\n")
  cat("\n")
  cat("index   ", y$first, ":", y$last, "\n", sep = "")
  cat(" date   (", as.character(y$date[y$first+1]), ", ", as.character(y$date[y$last+1]), "]\n", sep = "")
  cat("cases   ", sum(y$cases[y$first:y$last]), " of ", sum(y$cases), "\n", sep = "")
  cat("\n")
  cat("Fitted parameter estimates:\n")
  cat("\n")
  print(x$theta_fit)
  cat("\n")
  cat("Units:\n")
  cat("\n")
  print(uvec, quote = FALSE)
  cat("\n")
  cat("Negative log likelihood:", x$nll, "\n")
  invisible(x)
}

#' Print doubling times
#'
#' @description
#' A method for printing objects of class "doubling_time".
#'
#' @param x A "doubling_time" object.
#' @param ... Unused optional arguments.
#'
#' @return
#' `x` (invisibly).
#'
#' @seealso [compute_doubling_time()]
#' @export
print.doubling_time <- function(x, ...) {
  cat("Doubling times in days:\n")
  cat("\n")
  print(unclass(x))
  invisible(x)
}
