#' Extract fitted values
#'
#' Computes fitted values of nonlinear and dispersion model parameters.
#' The fitted value for a given fitting window is obtained by adding
#' (i) the population fitted value computed as a linear combination
#' of fixed effects coefficients and
#' (ii) all applicable random effects, with random effects coefficients
#' set equal to their conditional modes.
#'
#' @param object
#'   An \code{"\link{egf}"} object.
#' @param par
#'   A subset of \code{\link{get_par_names}(object, link = TRUE)}
#'   naming nonlinear and dispersion model parameters for which
#'   fitted values should be retrieved.
#' @param subset
#'   An expression to be evaluated in the combined model frame
#'   (see \code{\link{make_combined}}). Must evaluate to
#'   a \link{logical} vector indexing rows of the data frame,
#'   and thus fitting windows. Fitted values are retrieved
#'   only for indexed windows. The default (\code{\link{NULL}})
#'   is to consider all windows.
#' @param append
#'   An expression indicating variables in the combined model frame
#'   (see \code{\link{make_combined}}) to be included with the result.
#'   The default (\code{\link{NULL}}) is to append nothing.
#' @param link
#'   A \link{logical} flag. If \code{FALSE}, then fitted values
#'   are inverse link-transformed.
#' @param se
#'   A \link{logical} flag. If \code{link = TRUE} and \code{se = TRUE},
#'   then approximate (delta method) standard errors on fitted values
#'   are reported. Note that standard errors are required for subsequent
#'   use of \code{\link{confint.egf_fitted}}.
#' @param .subset
#'   A \link{logical} vector to be used (if non-\code{\link{NULL}})
#'   in place of the result of evaluating \code{subset}.
#' @param .append
#'   A \link{character} vector listing variable names to be used
#'   (if non-\code{\link{NULL}}) in place of the result of evaluating
#'   \code{append}.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' \code{coef.egf} is currently an alias for \code{fitted.egf}.
#'
#' See topic \code{\link{nse}} for details on nonstandard evaluation
#' of \code{subset} and \code{append}.
#'
#' @return
#' A \link[=data.frame]{data frame} inheriting from \link{class}
#' \code{"egf_fitted"}, with variables:
#' \item{par}{
#'   Nonlinear or dispersion model parameter,
#'   from \code{\link{get_par_names}(object, link = link)}.
#' }
#' \item{ts}{
#'   Time series, from \code{\link{levels}(object$endpoints$ts)}.
#' }
#' \item{window}{
#'   Fitting window, from \code{\link{levels}(object$endpoints$window)}.
#' }
#' \item{estimate}{
#'   Fitted value of parameter \code{par} in fitting window \code{window}.
#' }
#' \item{se}{
#'   (If \code{link = TRUE} and \code{se = TRUE}.)
#'   Approximate (delta method) standard error on \code{estimate}.
#' }
#'
#' @seealso \code{\link{confint.egf_fitted}}
#' @export
fitted.egf <- function(object,
                       par = get_par_names(object, link = TRUE),
                       subset = NULL,
                       append = NULL,
                       link = TRUE,
                       se = FALSE,
                       .subset = NULL,
                       .append = NULL,
                       ...) {
  stop_if_not_true_false(link)
  par_names <- get_par_names(object, link = TRUE)
  par <- unique(match.arg(par, par_names, several.ok = TRUE))
  combined <- make_combined(object)
  subset <- subset_to_index(substitute(subset), data = combined, enclos = parent.frame(),
                            .subset = .subset)
  append <- append_to_index(substitute(append), data = combined, enclos = parent.frame(),
                            .append = .append)

  if (link && se) {
    if (is.null(object$sdreport)) {
      warning(wrap(
        "Computing a Hessian matrix, which could take a while. ",
        "To avoid needless recomputation, do ",
        "`object$sdreport <- try(TMB::sdreport(object$tmb_out))` first."
      ))
      object$sdreport <- try(TMB::sdreport(object$tmb_out))
    }
    if (inherits(object$sdreport, "try-error")) {
      stop(wrap(
        "Unable to proceed because `TMB::sdreport(object$tmb_out)` ",
        "throws an error. Retry after diagnosing and refitting."
      ))
    }
    Y <- as.list(object$sdreport, what = "Estimate", report = TRUE)$Y_as_vector
    Y_se <- as.list(object$sdreport, what = "Std. Error", report = TRUE)$Y_as_vector
  } else {
    Y <- object$tmb_out$report(object$best)$Y_as_vector
  }

  ## `Y[i, j]` is the fitted value of nonlinear or dispersion parameter `j`
  ## (link scale) in fitting window `i`
  dim(Y) <- c(nrow(object$endpoints), length(par_names))
  colnames(Y) <- par_names
  Y <- Y[subset, par, drop = FALSE]

  d <- data.frame(
    par = rep(factor(par, levels = par_names), each = sum(subset)),
    object$endpoints[subset, c("ts", "window"), drop = FALSE],
    estimate = as.numeric(Y)
  )
  if (link && se) {
    dim(Y_se) <- dim(Y)
    colnames(Y_se) <- colnames(Y)
    Y_se <- Y_se[subset, par, drop = FALSE]
    d$se <- as.numeric(Y_se)
  }
  if (!link) {
    d$estimate <- mftapply(d$estimate, d$par,
      f = lapply(string_extract_link(levels(d$par)), match_link, inverse = TRUE)
    )
    levels(d$par) <- string_remove_link(levels(d$par))
  }
  d <- data.frame(
    d,
    combined[subset, append, drop = FALSE],
    row.names = NULL,
    check.names = FALSE
  )
  attr(d, "se") <- link && se
  class(d) <- c("egf_fitted", "data.frame")
  d
}

#' @rdname fitted.egf
#' @export
coef.egf <- fitted.egf

#' Confidence intervals on fitted values
#'
#' Computes confidence intervals on fitted values of nonlinear and dispersion
#' model parameters.
#'
#' @param object
#'   An \code{"\link[=fitted.egf]{egf_fitted}"} object.
#'   Must supply standard errors on fitted values.
#' @param parm
#'   Unused argument included for generic consistency.
#' @param level
#'   A number in the interval (0,1) indicating a confidence level.
#' @param link
#'   A \link{logical} flag. If \code{FALSE}, then confidence intervals
#'   on inverse link-transformed fitted values are returned.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' Confidence limits on fitted values (link scale) are computed
#' as \code{estimate + c(-1, 1) * sqrt(q) * se},
#' with \code{estimate} and \code{se} as in \code{object} and
#' \code{q = \link{qchisq}(level, df = 1)}.
#'
#' @return
#' If \code{link = TRUE}, then \code{object} but with variable
#' \code{se} replaced with variables \code{lower} and \code{upper}
#' supplying confidence limits on fitted values (link scale).
#'
#' Otherwise, the same object but with variables \code{estimate},
#' \code{lower}, and \code{upper} inverse link-transformed and
#' \code{\link{levels}(par)} modified accordingly.
#'
#' \code{level} is retained as an \link[=attributes]{attribute}.
#'
#' @export
confint.egf_fitted <- function(object, parm, level = 0.95, link = TRUE, ...) {
  stop_if_not(
    attr(object, "se"),
    m = wrap(
      "`object` must supply standard errors on fitted values. ",
      "Repeat `fitted()` with `link = TRUE` and `se = TRUE`."
    )
  )
  stop_if_not_number_in_interval(level, 0, 1, "()")
  stop_if_not_true_false(link)

  s <- c("par", "ts", "window", "estimate", "se")
  d <- data.frame(
    object[s[1:4]],
    do_wald(estimate = object$estimate, se = object$se, level = level),
    object[-match(s, names(object), 0L)],
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  attr(d, "level") <- level
  if (link) {
    return(d)
  }
  s_elu <- c("estimate", "lower", "upper")
  d[s_elu] <- mftapply(d[s_elu], d$par,
    f = lapply(string_extract_link(levels(d$par)), match_link, inverse = TRUE)
  )
  levels(d$par) <- string_remove_link(levels(d$par))
  d
}
