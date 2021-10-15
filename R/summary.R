#' Summarize a model object
#'
#' Summarizes fitted values of top level nonlinear model parameters
#' and gathers diagnostic information that can be used to quickly
#' assess convergence of the optimizer.
#'
#' @param object
#'   An \code{"\link{egf}"} object.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A list inheriting from class \code{"egf_summary"}, with elements:
#' \item{fitted}{
#'   A numeric matrix. Each column is the result of applying
#'   \code{\link{summary.default}} to a numeric vector listing
#'   the fitted values of a top level nonlinear model parameter.
#'   Fitted values are obtained via \code{\link[=fitted.egf]{fitted}}.
#' }
#' \item{convergence}{
#'   An integer code returned by the optimizer, with 0 indicating successful
#'   convergence within the specified absolute or relative tolerance.
#' }
#' \item{value, gradient}{
#'   Numeric vectors giving the value and gradient of the negative log
#'   likelihood function at the parameter vector returned by the optimizer.
#' }
#' \item{hessian}{
#'   A logical flag indicating whether the Hessian matrix of the negative log
#'   likelihood function at the parameter vector returned by the optimizer is
#'   positive definite.
#'   \code{NA} means that the matrix was not computed by \code{\link{egf}},
#'   either because \code{\link{egf}} was not called with \code{se = TRUE},
#'   or because an error was thrown during computation.
#'   In the first case, \code{object$sdreport} is \code{NULL}.
#'   In the second case, it is a \code{"\link[=try]{try-error}"} object
#'   preserving the error message.
#' }
#'
#' @export
#' @importFrom stats fitted aggregate
summary.egf <- function(object, ...) {
  safe_summary <- function(x) {
    res <- unclass(summary(x))
    if (anyNA(x)) {
      return(res)
    }
    c(res, `NA's` = 0)
  }
  fo <- fitted(object)
  sfo <- aggregate(fo["estimate"], fo["top"], safe_summary)
  sfo <- `colnames<-`(t(sfo[[2L]]), as.character(sfo[[1L]]))
  if (sum(sfo["NA's", ] == 0)) {
    sfo <- sfo[-match("NA's", rownames(sfo), 0L), , drop = FALSE]
  }
  res <- list(
    fitted = sfo,
    convergence = object$optimizer_out$convergence,
    value = object$value,
    gradient = object$gradient,
    hessian = object$hessian
  )
  class(res) <- "egf_summary"
  res
}

#' @export
print.egf_summary <- function(x, width = 0.9 * getOption("width"), indent = 2L, ...) {
  heading("Fitted values", width = width)
  cat("\n")
  writeLines(paste0(strrep(" ", indent), capture.output(print(x$fitted))))
  cat("\n")

  heading("Negative log likelihood", width = width)
  cat("\n")
  c1 <- c(
    "convergence",
    "value",
    "range(abs(gradient))",
    "pos. def. Hessian"
  )
  c2 <- c(
    sprintf("%d", x$convergence),
    sprintf("%.6e", x$value),
    paste0(sprintf("%.6e", range(abs(x$gradient))), collapse = " "),
    as.character(x$hessian)
  )
  writeLines(paste0(strrep(" ", indent), align(c1, c2, justify = "l")))

  invisible(x)
}
