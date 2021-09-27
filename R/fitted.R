#' Compute fitted values
#'
#' Computes fitted values of top level nonlinear model parameters.
#' The fitted value of a given parameter for a given fitting window
#' is obtained by adding
#' (i) the population fitted value computed as a linear combination
#' of fixed effects coefficients and
#' (ii) all applicable random effects, with random effects coefficients
#' set equal to their conditional modes.
#'
#' @param object
#'   An \code{"\link{egf}"} object.
#' @param top
#'   A subset of \code{\link{egf_get_names_top}(object, link = TRUE)}
#'   naming top level nonlinear model parameters for which fitted
#'   values should be retrieved.
#' @param subset
#'   An expression to be evaluated in the combined model frame
#'   (see \code{\link{egf_combine_frames}}). It must evaluate
#'   to a valid index vector for the rows of the data frame
#'   (see \code{\link{[.data.frame}}), and thus fitting windows.
#'   Fitted values are retrieved only for indexed windows.
#'   The default (\code{\link{NULL}}) is to consider all windows.
#' @param append
#'   An expression indicating variables in the combined model frame
#'   (see \code{\link{egf_combine_frames}}) to be included with the
#'   result. The default (\code{\link{NULL}}) is to append nothing.
#' @param link
#'   A \link{logical} flag. If \code{FALSE}, then fitted values are
#'   inverse link-transformed.
#' @param se
#'   A \link{logical} flag. If \code{link = TRUE} and \code{se = TRUE},
#'   then approximate delta method standard errors on fitted values
#'   are reported. Standard errors are required for subsequent use of
#'   \code{\link{confint.egf_fitted}}.
#' @param .subset,.append
#'   Index vectors to be used (if non-\code{\link{NULL}}) in place of
#'   the result of evaluating \code{subset} and \code{append}.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' \code{coef.egf} is currently an alias for \code{fitted.egf}.
#'
#' See topic \code{\link{egf_eval}} for details on nonstandard evaluation
#' of \code{subset} and \code{append}.
#'
#' @return
#' A \link[=data.frame]{data frame} inheriting from \link{class}
#' \code{"egf_fitted"}, with variables:
#' \item{top}{
#'   Top level nonlinear model parameter,
#'   from \code{\link{egf_get_names_top}(object, link = link)}.
#' }
#' \item{ts}{
#'   Time series, from \code{\link{levels}(object$frame_windows$ts)}.
#' }
#' \item{window}{
#'   Fitting window, from \code{\link{levels}(object$frame_windows$window)}.
#' }
#' \item{estimate}{
#'   Fitted value of parameter \code{top} in fitting window \code{window}.
#' }
#' \item{se}{
#'   (If \code{link = TRUE} and \code{se = TRUE}.)
#'   Approximate delta method standard error on \code{estimate}.
#' }
#'
#' @examples
#' example("egf", package = "epigrowthfit", local = TRUE, echo = FALSE)
#' exdata <- system.file("exdata", package = "epigrowthfit", mustWork = TRUE)
#' object <- readRDS(file.path(exdata, "egf.rds"))
#'
#' path_to_cache <- file.path(exdata, "egf_fitted.rds")
#' if (file.exists(path_to_cache)) {
#'   zz <- readRDS(path_to_cache)
#' } else {
#'   zz <- fitted(object, se = TRUE)
#'   saveRDS(zz, file = path_to_cache)
#' }
#' str(zz)
#'
#' @family coefficient extractors
#' @seealso \code{\link{confint.egf_fitted}}
#' @export
fitted.egf <- function(object,
                       top = egf_get_names_top(object, link = TRUE),
                       subset = NULL,
                       append = NULL,
                       link = TRUE,
                       se = FALSE,
                       .subset = NULL,
                       .append = NULL,
                       ...) {
  stop_if_not_true_false(link)
  stop_if_not_true_false(se)
  names_top <- egf_get_names_top(object, link = TRUE)
  top <- unique(match.arg(top, names_top, several.ok = TRUE))

  combined <- egf_combine_frames(object)
  subset <- if (is.null(.subset)) substitute(subset) else .subset
  subset <- egf_eval_subset(subset, combined, parent.frame())
  append <- if (is.null(.append)) substitute(append) else .append
  append <- egf_eval_append(append, combined, baseenv())

  if (link && se) {
    if (is.null(object$sdreport)) {
      if (egf_has_random(object)) {
        warning(wrap(
          "Computing a Hessian matrix for a model with random effects, ",
          "which might take a while. To avoid needless recomputation, ",
          "do 'object$sdreport <- try(TMB::sdreport(object$tmb_out))' ",
          "before trying 'fitted'."
        ))
      }
      object$sdreport <- try(TMB::sdreport(object$tmb_out), silent = TRUE)
    }
    if (inherits(object$sdreport, "try-error")) {
      stop(wrap(
        "Unable to proceed because 'TMB::sdreport(object$tmb_out)' ",
        "throws the following error:\n\n",
        conditionMessage(attr(object$sdreport, "condition")), "\n\n",
        "Retry after diagnosing and refitting."
      ))
    }
    ssdr <- summary(object$sdreport, select = "report")
    index <- rownames(ssdr) == "Y"
    Y <- ssdr[index, "Estimate"]
    Y_se <- ssdr[index, "Std. Error"]
    dim(Y) <- dim(Y_se) <- object$tmb_out$env$ADreportDims$Y
  } else {
    Y <- object$tmb_out$report(object$best)$Y
  }

  ## 'Y[i, j]' is the fitted value of top level nonlinear model parameter 'j'
  ## (link scale) in fitting window 'i'
  colnames(Y) <- names_top
  Y <- Y[subset, top, drop = FALSE]

  res <- data.frame(
    top = rep(factor(top, levels = names_top), each = length(subset)),
    object$frame_windows[subset, c("ts", "window"), drop = FALSE],
    estimate = as.numeric(Y),
    row.names = NULL,
    check.names = FALSE
  )
  if (link && se) {
    colnames(Y_se) <- names_top
    Y_se <- Y_se[subset, top, drop = FALSE]
    res$se <- as.numeric(Y_se)
  }
  if (!link) {
    res$estimate <- in_place_ragged_apply(res$estimate, res$top,
      f = lapply(egf_link_extract(levels(res$top)), egf_link_match, inverse = TRUE)
    )
    levels(res$top) <- egf_link_remove(levels(res$top))
  }
  res <- data.frame(
    res,
    combined[subset, append, drop = FALSE],
    row.names = NULL,
    check.names = FALSE
  )
  attr(res, "se") <- link && se
  class(res) <- c("egf_fitted", "data.frame")
  res
}

#' @rdname fitted.egf
#' @export
coef.egf <- fitted.egf

#' Confidence intervals on fitted values
#'
#' Computes confidence intervals on fitted values of top level nonlinear model
#' parameters.
#'
#' @param object
#'   An \code{"\link[=fitted.egf]{egf_fitted}"} object.
#'   Must supply link scale fitted values and corresponding standard errors.
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
#' as \code{estimate[i] + c(-1, 1) * sqrt(q) * se[i]},
#' with \code{estimate} and \code{se} as in \code{object} and
#' \code{q = \link{qchisq}(level, df = 1)}.
#'
#' @return
#' If \code{link = TRUE}, then \code{object} but with variable
#' \code{se} replaced with variables \code{lower} and \code{upper}
#' supplying confidence limits on fitted values (link scale).
#'
#' Otherwise, the same result but with variables \code{estimate},
#' \code{lower}, and \code{upper} inverse link-transformed and
#' the \link{levels} of variable \code{top} modified accordingly.
#'
#' \code{level} is retained as an \link[=attributes]{attribute}.
#'
#' @examples
#' example("fitted.egf", package = "epigrowthfit", local = TRUE, echo = FALSE)
#' object <- readRDS(system.file("exdata", "egf_fitted.rds",
#'                               package = "epigrowthfit", mustWork = TRUE))
#'
#' confint(object, link = TRUE)
#' confint(object, link = FALSE)
#'
#' @export
confint.egf_fitted <- function(object, parm, level = 0.95, link = TRUE, ...) {
  if (!isTRUE(attr(object, "se"))) {
    stop(wrap(
      "'object' must supply link scale fitted values ",
      "and corresponding standard errors. ",
      "Retry with 'object = fitted(., link = TRUE, se = TRUE)'."
    ))
  }
  stop_if_not_number_in_interval(level, 0, 1, "()")
  stop_if_not_true_false(link)

  s <- c("top", "ts", "window", "estimate", "se")
  res <- data.frame(
    object[s[1:4]],
    do_wald(estimate = object$estimate, se = object$se, level = level),
    object[-match(s, names(object), 0L)],
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!link) {
    elu <- c("estimate", "lower", "upper")
    res[elu] <- in_place_ragged_apply(res[elu], res$top,
      f = lapply(egf_link_extract(levels(res$top)), egf_link_match, inverse = TRUE)
    )
    levels(res$top) <- egf_link_remove(levels(res$top))
  }
  attr(res, "level") <- level
  res
}
