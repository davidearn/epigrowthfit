#' Print models of epidemic growth
#'
#' @description
#' Methods for printing objects of class "egf_init" or "egf".
#'
#' @param x
#'   An "egf_init" or "egf" object.
#' @param ...
#'   Unused optional arguments.
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
  ### SET UP ###########################################################

  cstr <- switch(x$curve,
    exponential = "an exponential model",
    logistic    = "a logistic model",
    richards    = "a Richards model"
  )
  bstr <- if (x$include_baseline) {
    "a linear baseline and "
  } else {
    ""
  }
  dstr <- switch(x$distr,
    poisson = "Poisson-distributed",
    nbinom  = "negative binomial"
  )
  mstr <- sprintf(
    fmt = paste0(
      "Pass this \"egf_init\" object to `egf()` to fit %s\n",
      "with %s%s observations."
    ),
    cstr, bstr, dstr
  )

  i1 <- x$first
  i2 <- x$last
  d1 <- x$date[i1]
  d2 <- x$date[i2]
  c1 <- sum(x$cases[(i1+1):i2], na.rm = TRUE)
  c2 <- sum(x$cases[-1], na.rm = TRUE)
  nna <- sum(is.na(x$cases[(i1+1):i2]))
  wstr <- sprintf(
    fmt = paste0(
      "index   %d:%d\n",
      " date   %s to %s (%d days)\n",
      "cases   %d of %d\n",
      "   NA   %d"
    ),
    i1, i2,
    as.character(d1), as.character(d2), days(d2, since = d1),
    c1, c2,
    nna
  )

  pvec <- x$theta_init
  uvec <- c(r = "per day", thalf = "days", b = "per day")
  uvec <- uvec[names(uvec) %in% names(pvec)]


  ### PRINT ############################################################

  cat(mstr, "\n")
  cat("\n")
  cat("Fitting window:\n")
  cat("\n")
  cat(wstr, "\n")
  cat("\n")
  cat("Initial parameter estimates:\n")
  cat("\n")
  print(pvec)
  cat("\n")
  cat("Units:\n")
  cat("\n")
  print(uvec, quote = FALSE)
  invisible(x)
}

#' @rdname print.egf
#' @export
print.egf <- function(x, ...) {
  ### SET UP ###########################################################

  y <- x$init
  cstr <- switch(y$curve,
    exponential = "an exponential model",
    logistic    = "a logistic model",
    richards    = "a Richards model"
  )
  bstr <- if (y$include_baseline) {
    "a linear baseline and "
  } else {
    ""
  }
  dstr <- switch(y$distr,
    poisson = "Poisson-distributed",
    nbinom  = "negative binomial"
  )
  mstr <- sprintf(
    fmt = paste0(
      "This \"egf\" object fits %s\n",
      "with %s%s observations."
    ),
    cstr, bstr, dstr
  )

  i1 <- y$first
  i2 <- y$last
  d1 <- y$date[i1]
  d2 <- y$date[i2]
  c1 <- sum(y$cases[(i1+1):i2], na.rm = TRUE)
  c2 <- sum(y$cases[-1], na.rm = TRUE)
  nna <- sum(is.na(y$cases[(i1+1):i2]))
  wstr <- sprintf(
    fmt = paste0(
      "index   %d:%d\n",
      " date   %s to %s (%d days)\n",
      "cases   %d of %d\n",
      "   NA   %d"
    ),
    i1, i2,
    as.character(d1), as.character(d2), days(d2, since = d1),
    c1, c2,
    nna
  )

  pvec <- x$theta_fit
  uvec <- c(r = "per day", thalf = "days", b = "per day")
  uvec <- uvec[names(uvec) %in% names(pvec)]


  ### PRINT ############################################################

  cat(mstr, "\n")
  cat("\n")
  cat("Fitting window:\n")
  cat("\n")
  cat(wstr, "\n")
  cat("\n")
  cat("Fitted parameter estimates:\n")
  cat("\n")
  print(pvec)
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
#' @param x
#'   A "doubling_time" object.
#' @param ...
#'   Unused optional arguments.
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
