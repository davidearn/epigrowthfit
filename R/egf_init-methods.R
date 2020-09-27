#' Methods for class "epigrowthfit_init"
#'
#' @description
#' Methods for printing and plotting "egf_init" objects
#' returned by [egf_init()].
#'
#' @param x An "epigrowthfit_init" object.
#' @param ... Unused optional arguments.
#'
#' @name egf_init-methods
NULL

#' @rdname egf_init-methods
#' @export
#' @import graphics
print.egf_init <- function(x, ...) {
  if (!inherits(x, "egf_init")) {
    stop("`x` must be an \"egf_init\" object.")
  }
  model <- switch(x$arg_list$model,
    exponential = "an exponential",
    logistic    = "a logistic",
    richards    = "a Richards"
  )
  distribution <- switch(x$arg_list$distribution,
    poisson = "Poisson-distributed",
    nbinom  = "negative binomial"
  )
  first <- x$first
  last <- x$last
  time <- x$arg_list$time
  theta0 <- x$theta0
  cat("Pass this \"egf_init\" object to `egf()` to fit\n")
  cat(model, "model assuming", distribution, "observations.\n")
  cat("\n")
  cat("Fitting window:\n")
  cat("\n")
  cat("index   ", first, ":", last, "\n", sep = "")
  cat(" time   [", time[first], ", ", time[last], "]\n", sep = "")
  cat("\n")
  cat("Initial parameter estimates:\n")
  cat("\n")
  print(unlist(theta0))
  invisible(x)
}

#' @rdname egf_init-methods
#' @export
#' @import graphics
plot.egf_init <- function(x, ...) {
  if (!inherits(x, "egf_init")) {
    stop("`x` must be an \"egf_init\" object.")
  }
  op <- par(mar = c(5, 4, 4, 8) + 0.1)
  time <- x$arg_list$time
  cases <- x$arg_list$cases
  first <- x$first
  last <- x$last
  theta0 <- x$theta0
  plot(log(cases + 0.1) ~ time, data = data.frame(time, cases), las = 1)
  abline(v = time[c(first, last)], lty = 2)
  axis(side = 3, at = time[c(first, last)], labels = c(first, last),
       tick = FALSE, mgp = c(3, 0.1, 0))
  mtext("index", side = 3, line = 2)
  ## Hack to get alignment at "=" and at "e"
  ## when listing initial parameter estimates
  str1 <- paste0(names(theta0), " = ")
  mat <- matrix(unlist(strsplit(sprintf("%0.3e", unlist(theta0)), "e")), nrow = 2)
  str2 <- paste0(mat[1, ], "e")
  str3 <- paste0(mat[2, ])
  xx1 <- par("usr")[2] + 0.02 * diff(par("usr")[1:2]) + max(strwidth(str1))
  xx2 <- xx1 + max(strwidth(str2))
  yy <- mean(par("usr")[3:4])
  text(xx1, yy, paste(str1, collapse = "\n"), adj = c(1, 0.5), xpd = NA)
  text(xx2, yy, paste(str2, collapse = "\n"), adj = c(1, 0.5), xpd = NA)
  text(xx2, yy, paste(str3, collapse = "\n"), adj = c(0, 0.5), xpd = NA)
  par(op)
  invisible(NULL)
}

