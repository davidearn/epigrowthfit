#' Extract fitted values
#'
#' Retrieves fitted values of top level nonlinear model parameters.
#' The fitted value of a given parameter for a given fitting window
#' is obtained by adding
#' (i) the population fitted value computed as a linear combination
#' of fixed effects coefficients and
#' (ii) all applicable random effects, with random effects set equal
#' to their conditional modes.
#'
#' @param object
#'   An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.
#' @param top
#'   A subset of \code{\link{egf_get_names_top}(object, link = TRUE)}
#'   naming top level nonlinear model parameters for which fitted
#'   values should be retrieved.
#' @param link
#'   A logical flag.
#'   If \code{FALSE}, then fitted values are inverse link-transformed.
#' @param se
#'   A logical flag. If \code{se = TRUE} and \code{link = TRUE},
#'   then approximate delta method standard errors on fitted values
#'   are reported.
#'   Standard errors are required for subsequent use
#'   of \code{\link{confint.egf_fitted}}.
#'   Setting \code{se = TRUE} and \code{link = FALSE} is an error,
#'   as standard errors are not available for inverse link-transformed
#'   fitted values.
#'   Setting \code{se = TRUE} when \code{object} inherits from class
#'   \code{"egf_no_fit"} is likewise an error.
#' @param subset
#'   An expression to be evaluated in
#'   \code{\link[=model.frame.egf]{model.frame}(object, "combined")}.
#'   It must evaluate to a valid index vector for the rows of
#'   the data frame and, in turn, fitting windows.
#'   Fitted values are retrieved only for indexed windows.
#'   The default (\code{NULL}) is to consider all windows.
#' @param append
#'   An expression indicating variables in
#'   \code{\link[=model.frame.egf]{model.frame}(object, "combined")}
#'   to be included with the result.
#'   The default (\code{NULL}) is to append nothing.
#' @param .subset,.append
#'   Index vectors to be used (if non-\code{NULL}) in place of
#'   the result of evaluating \code{subset} and \code{append}.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' See topic \code{\link{egf_eval}} for details on nonstandard evaluation
#' of \code{subset} and \code{append}.
#'
#' @return
#' A data frame inheriting from class \code{"egf_fitted"}, with variables:
#' \item{top}{
#'   Top level nonlinear model parameter,
#'   from \code{\link{egf_get_names_top}(object, link = link)}.
#' }
#' \item{ts}{
#'   Time series, from
#'   \code{levels(\link[=model.frame.egf]{model.frame}(object)$ts)}.
#' }
#' \item{window}{
#'   Fitting window, from
#'   \code{levels(\link[=model.frame.egf]{model.frame}(object)$window)}.
#' }
#' \item{estimate}{
#'   Fitted value of parameter \code{top} in fitting window \code{window}.
#'   If \code{object} inherits from class \code{"egf"},
#'   then this is conditioned on the estimated mixed effects model.
#'   If \code{object} inherits from class \code{"egf_no_fit"},
#'   then this is conditioned on the \emph{initial} values
#'   of mixed effects model parameters.
#' }
#' \item{se}{
#'   Approximate delta method standard error on \code{estimate}.
#'   Absent except for calls matching \code{fitted(link = TRUE, se = TRUE)}.
#' }
#'
#' @examples
#' object <- egf_cache("egf-1.rds")
#' zz <- egf_cache("fitted-egf-1.rds", fitted(object, se = TRUE))
#' str(zz)
#'
#' @family extractors
#' @seealso \code{\link{confint.egf_fitted}}
#' @export
#' @importFrom stats model.frame
#' @importFrom TMB sdreport
fitted.egf <- function(object,
                       top = egf_get_names_top(object, link = TRUE),
                       link = TRUE,
                       se = FALSE,
                       subset = NULL,
                       append = NULL,
                       .subset = NULL,
                       .append = NULL,
                       ...) {
  stopifnot(
    is_true_or_false(link),
    is_true_or_false(se)
  )
  if (se && !link) {
    stop(wrap(
      "Standard errors are not available for inverse link-transformed ",
      "fitted values."
    ))
  }

  names_top <- egf_get_names_top(object, link = TRUE)
  top <- unique(match.arg(top, names_top, several.ok = TRUE))

  frame_windows <- model.frame(object, "windows")
  frame_combined <- model.frame(object, "combined")
  subset <- if (is.null(.subset)) substitute(subset) else .subset
  subset <- egf_eval_subset(subset, frame_combined, parent.frame())
  append <- if (is.null(.append)) substitute(append) else .append
  append <- egf_eval_append(append, frame_combined, baseenv())

  if (se) {
    sdr <- egf_get_sdreport(object)
    ssdr <- summary(sdr, select = "report")
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
    frame_windows[subset, c("ts", "window"), drop = FALSE],
    estimate = as.numeric(Y),
    row.names = NULL,
    check.names = FALSE
  )
  if (se) {
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
    frame_combined[subset, append, drop = FALSE],
    row.names = NULL,
    check.names = FALSE
  )
  attr(res, "se") <- se
  class(res) <- c("egf_fitted", oldClass(res))
  res
}

#' @rdname fitted.egf
#' @export
fitted.egf_no_fit <- function(object,
                              top = egf_get_names_top(object, link = TRUE),
                              link = TRUE,
                              se = FALSE,
                              subset = NULL,
                              append = NULL,
                              .subset = NULL,
                              .append = NULL,
                              ...) {
  if (se) {
    stop(wrap(
      "Standard errors cannot be computed until the model is estimated. ",
      "Retry after doing, e.g., 'object <- update(object, se = TRUE, fit = TRUE, ...)'."
    ))
  }

  ## Passing arguments to method for class "egf" without evaluating
  ## 'subset' or 'append' requires minor acrobatics
  call <- match.call(expand.dots = FALSE)
  call[[1L]] <- quote(fitted.egf)
  call$... <- NULL
  nms <- names(call)
  i <- match(nms[-1L], c("subset", "append"), 0L) == 0L
  call[-1L][i] <- lapply(nms[-1L][i], as.name)

  object$best <- object$init
  eval(call)
}

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
#' object <- egf_cache("fitted-egf-1.rds")
#' zz <- confint(object)
#' str(zz)
#'
#' @export
confint.egf_fitted <- function(object, parm, level = 0.95, link = TRUE, ...) {
  if (!isTRUE(attr(object, "se"))) {
    stop(wrap(
      "'object' must supply link scale fitted values ",
      "and corresponding standard errors. ",
      "Retry with 'object = fitted(<\"egf\" object>, link = TRUE, se = TRUE)'."
    ))
  }
  stopifnot(
    is_number_in_interval(level, 0, 1, "()"),
    is_true_or_false(link)
  )

  s <- c("top", "ts", "window", "estimate", "se")
  res <- data.frame(
    object[s[1:4]],
    wald(estimate = object$estimate, se = object$se, level = level),
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
