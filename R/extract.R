#' Extract coefficients and covariance parameters
#'
#' Extracts bottom level parameter vectors \code{beta}, \code{theta},
#' and \code{b}. These store
#' linear fixed effects coefficients,
#' random effect covariance parameters, and
#' linear random effects coefficients, respectively.
#'
#' @param object
#'   An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.
#' @param full
#'   A logical flag. If \code{FALSE}, then parameter vectors are returned
#'   in the condensed format used by \pkg{TMB}, which excludes mapped elements.
#'   (See \code{\link{egf}} argument \code{map}.)
#' @param list
#'   A logical flag.
#'   If \code{FALSE}, then the parameter vectors are concatenated.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' If \code{object} inherits from class \code{"egf_no_fit"},
#' then the result lists initial parameter values to be used
#' in the first likelihood evaluation.
#'
#' @return
#' If \code{list = FALSE}, then a numeric vector \code{x = c(beta, theta, b)},
#' with \code{names(x)} providing the grouping.
#'
#' If \code{list = TRUE}, then a named list \code{list(beta, theta, b)}.
#' Each parameter vector has an attribute \code{map}, an index vector such
#' that the full vector \code{y} and condensed vector \code{x} are related
#' by \code{y = x[map]}, with the exception that \code{map[i]} is \code{NA}
#' if \code{y[i]} is mapped to an initial value
#' (hence \code{x[map[i]]} is equal to \code{NA}, not \code{y[i]}).
#'
#' In both cases, the result inherits from class \code{"egf_coef"}.
#'
#' @family extractors
#' @export
#' @importFrom stats coef
coef.egf <- function(object, full = FALSE, list = FALSE, ...) {
  stop_if_not_true_false(full)
  stop_if_not_true_false(list)

  if (full) {
    res <- egf_expand_par(object$tmb_out, par = object$best)
  } else {
    res <- object$best
    attr(res, "lengths") <- c(table(factor(names(res), levels = c("beta", "theta", "b"))))
  }
  if (list) {
    names(res) <- NULL
    len <- attr(res, "lengths")
    res <- split(res, rep.int(gl(length(len), 1L, labels = names(len)), len))
    map <- lapply(object$tmb_out$env$parList(), attr, "map")
    for (s in names(res)) {
      if (length(res[[s]]) == 0L || is.null(map[[s]])) {
        attr(res[[s]], "map") <- seq_along(res[[s]])
      } else {
        map[[s]] <- map[[s]] + 1L
        map[[s]][map[[s]] == 0L] <- NA
        attr(res[[s]], "map") <- map[[s]]
      }
    }
  }
  attr(res, "full") <- full
  class(res) <- "egf_coef"
  res
}

#' @rdname coef.egf
#' @export
#' @importFrom stats coef
coef.egf_no_fit <- function(object, full = FALSE, list = FALSE, ...) {
  object$best <- object$init
  coef.egf(object, full = full, list = list, ...)
}

#' @export
print.egf_coef <- function(x, ...) {
  y <- x
  if (is.list(x)) {
    attr(x, "lengths") <- NULL
    for (i in seq_along(x)) {
      attr(x[[i]], "map") <- NULL
    }
  }
  attr(x, "full") <- NULL
  class(x) <- NULL
  NextMethod("print")
  invisible(y)
}

#' @method as.data.frame egf_coef
#' @export
as.data.frame.egf_coef <- as.data.frame.vector

#' Extract fixed effect coefficients
#'
#' Retrieve the coefficients of the fixed effects component
#' of a mixed effects model.
#'
#' @param object
#'   An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A data frame with one row per coefficient and variables:
#' \item{bottom}{
#'   Name of a bottom level mixed effects model parameter,
#'   in this case a fixed effects coefficient;
#'   this is a string of the form \code{"beta[\%d]"}.
#' }
#' \item{top}{
#'   Name of the top level nonlinear model parameter whose fitted value
#'   is a function of \code{bottom},
#'   from \code{\link{egf_get_names_top}(object, link = TRUE)}.
#' }
#' \item{term}{
#'   Term from the fixed effects component of the mixed effects model formula
#'   for parameter \code{top}.
#' }
#' \item{colname}{
#'   Column name in the fixed effects design matrix
#'   \code{\link[=model.matrix.egf]{model.matrix}(object, "fixed")}.
#' }
#' \item{estimate}{
#'   Coefficient estimate, from segment \code{beta} of
#'   \code{\link[=coef.egf]{coef}(object, full = TRUE)}.
#' }
#'
#' @family extractors
#' @aliases fixef
#' @export
#' @export fixef
#' @importFrom nlme fixef
fixef.egf <- function(object, ...) {
  par <- coef(object, full = TRUE)
  data.frame(
    object$effects$X[c("bottom", "top", "term", "colname")],
    estimate = par[names(par) == "beta"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' @rdname fixef.egf
#' @export
#' @importFrom nlme fixef
fixef.egf_no_fit <- fixef.egf

#' Extract random effect conditional modes
#'
#' Retrieve the coefficients of the random effects component
#' of a mixed effects model
#' (specifically, their modes conditional on data and parameter estimates).
#'
#' @param object
#'   An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.
#' @param build_cov
#'   A logical flag.
#'   If \code{TRUE}, then random effect covariance matrices are constructed
#'   from segment \code{theta} of
#'   \code{\link[=coef.egf]{coef}(object, full = TRUE)}
#'   and preserved as an attribute of the result.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A data frame with one row per coefficient and variables:
#' \item{cov}{
#'   Name of a covariance matrix.
#'   This is the interaction of \code{term} and \code{group},
#'   but with levels named \code{"Sigma[\%d]"}.
#' }
#' \item{vec}{
#'   Name of a random vector.
#'   This is the interaction of \code{term}, \code{group}, and \code{level},
#'   but with levels named \code{"u[\%d]"}.
#' }
#' \item{bottom}{
#'   Name of a bottom level mixed effects model parameter,
#'   in this case a random effects coefficient;
#'   this is a string of the form \code{"b[\%d]"}.
#' }
#' \item{top}{
#'   Name of the top level nonlinear model parameter whose fitted value
#'   is a function of \code{bottom}.
#' }
#' \item{term, group}{
#'   Random effects term from mixed effects model formula
#'   for parameter \code{top}. \code{term} and \code{group} give
#'   the left and right hand sides of the \code{`|`} operator.
#' }
#' \item{level}{
#'   Level of factor or interaction indicated by \code{group}.
#' }
#' \item{colname}{
#'   Column name in the random effects design matrix
#'   \code{\link[=model.matrix.egf]{model.matrix}(object, "random")}.
#' }
#' \item{mode}{
#'   Random effect conditional mode (unit variance scale),
#'   from segment \code{b} of
#'   \code{\link[=coef.egf]{coef}(object, full = TRUE)}.
#' }
#' If \code{build_cov = TRUE}, then the result has attribute \code{Sigma},
#' a list of covariance matrices corresponding to the levels of variable
#' \code{cov}.
#'
#' @family extractors
#' @aliases ranef
#' @export
#' @export ranef
#' @importFrom nlme ranef
ranef.egf <- function(object, build_cov = FALSE, ...) {
  stopifnot(egf_has_random(object))
  stop_if_not_true_false(build_cov)

  par <- coef(object, full = TRUE)
  res <- data.frame(
    object$effects$Z[c("cov", "vec", "bottom", "top", "term", "group", "level", "colname")],
    mode = par[names(par) == "b"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  if (!build_cov) {
    return(res)
  }

  d1 <- object$tmb_out$env$data$block_rows
  d2 <- object$tmb_out$env$data$block_cols
  p <- as.integer(choose(d1 + 1L, 2L))
  theta <- split(par[names(par) == "theta"], rep.int(seq_along(p), p))
  Sigma <- lapply(theta, theta2cov)
  names(Sigma) <- levels(object$effects$Z$cov)
  tt <- table(object$effects$Z$top, object$effects$Z$cov)
  for (i in seq_along(Sigma)) {
    dimnames(Sigma[[i]])[1:2] <- list(rownames(tt)[tt[, i] > 0L])
  }
  attr(res, "Sigma") <- Sigma
  res
}

#' @rdname ranef.egf
#' @export
#' @importFrom nlme ranef
ranef.egf_no_fit <- ranef.egf

#' Extract model covariance matrix
#'
#' Extracts (or, if necessary, computes) the covariance matrix
#' of bottom level parameters \code{beta} and \code{theta},
#' corresponding to the output of
#' \code{\link[=coef.egf]{coef}(object, full = FALSE)}.
#'
#' @param object
#'   An \code{"\link{egf}"} object.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' If \code{object} was constructed by a call \code{\link{egf}(se = TRUE)},
#' then the full covariance matrix has already been computed and is preserved
#' in \code{object}. \code{vcov} reuses this matrix to avoid recomputation.
#'
#' If the returned matrix is not finite and positive definite,
#' then the fit specified by \code{object} should be investigated,
#' as the optimizer that produced the fit may have failed to converge.
#' See also \code{\link{egf_has_converged}}.
#'
#' @return
#' A symmetric matrix.
#'
#' @family extractors
#' @export
#' @importFrom TMB sdreport
vcov.egf <- function(object, ...) {
  egf_get_sdreport(object)$cov.fixed
}

#' Extract model calls
#'
#' Extract the \link{call} to \code{\link{egf}} that produced the given
#' model object. These methods exist mainly to enable compatibility with
#' the default method of \code{\link{update}}.
#'
#' @param x
#'   An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A \link{call} to \code{\link{egf}}.
#'
#' @family extractors
#' @export
#' @importFrom stats getCall
getCall.egf <- function(x, ...) {
  call <- NextMethod("getCall")
  call[[1L]] <- quote(egf)
  call
}

#' @rdname getCall.egf
#' @export
#' @importFrom stats getCall
getCall.egf_no_fit <- getCall.egf

#' @export
#' @importFrom stats getCall
getCall.egf_model_simulate <- function(x, ...) {
  call <- NextMethod("getCall")
  call[[1L]] <- quote(simulate)
  call
}

#' Extract model frames
#'
#' Extract data frames, including \link[=model.frame]{model frame}s,
#' from a model object.
#'
#' @param formula
#'   An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.
#' @param which
#'   A character string controlling what is returned:
#'   \code{"ts"}, disease incidence time series;
#'   \code{"windows"}, fitting window endpoints;
#'   \code{"parameters"}, the mixed effects model frame
#'   corresponding to parameter \code{top};
#'   \code{"append"}, variables preserved in \code{formula}
#'   via \code{\link{egf}} argument \code{append}; or
#'   \code{"combined"}, the result of concatenating
#'   (in the sense of \code{\link{cbind}})
#'   the results of \code{"parameters"} and \code{"append"},
#'   then deleting any duplicated variables.
#' @param full
#'   A logical flag.
#'   If \code{TRUE}, then complete time series are returned.
#'   Otherwise, only observations belonging to fitting windows are returned.
#'   Ignored if \code{which != "ts"}.
#' @param top
#'   A character string specifying a top level nonlinear model parameter.
#'   Ignored if \code{which != "parameters"}.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' Only \code{which = "parameters"} results in a "proper" model frame,
#' i.e., a data frame with a \code{\link[=terms.object]{terms}} attribute.
#'
#' @return
#' A data frame.
#'
#' @family extractors
#' @method model.frame egf
#' @export
#' @importFrom stats model.frame
model.frame.egf <- function(formula,
                            which = c("ts", "windows", "parameters", "append", "combined"),
                            full = FALSE,
                            top = egf_get_names_top(formula, link = TRUE),
                            ...) {
  which <- match.arg(which)
  if (which == "combined") {
    res <- do.call(cbind, unname(formula$frame$parameters))
    res <- cbind(res, formula$frame$append)
    return(res[!duplicated(names(res))])
  }
  res <- formula$frame[[which]]
  if (which == "ts") {
    stop_if_not_true_false(full)
    if (full) {
      return(res)
    } else {
      return(res[!is.na(res$window), , drop = FALSE])
    }
  }
  if (which == "parameters") {
    top <- match.arg(top)
    return(res[[top]])
  }
  res
}

#' @rdname model.frame.egf
#' @method model.frame egf_no_fit
#' @export
#' @importFrom stats model.frame
model.frame.egf_no_fit <- model.frame.egf

#' Extract design matrices
#'
#' Extract fixed and random effects \link[=model.matrix]{design} matrices
#' from a model object.
#'
#' @param object
#'   An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.
#' @param which
#'   A character string controlling what is returned:
#'   \code{"fixed"}, the fixed effects design matrix corresponding
#'   to parameter \code{top}; or
#'   \code{"random"}, the random effects design matrix corresponding
#'   to parameter \code{top} and term \code{random}.
#' @param top
#'   A character string specifying a top level nonlinear model parameter.
#' @param random
#'   A random effects term, i.e., an expression of the form \code{(tt | g)}.
#'   Note that parentheses are required. Unused if \code{which = "fixed"}.
#' @param ...
#'   Unused optional arguments.
#'
#' @details
#' \code{model.matrix(which = "fixed", top = NULL)} returns the result of
#' combining (in the sense of \code{\link{cbind}}) all parameter-specific
#' fixed effects design matrices.
#'
#' \code{model.matrix(which = "random", top = <string>, random = NULL)}
#' returns the result of combining (in the sense of \code{\link{cbind}})
#' all term-specific random effects design matrices associated with
#' parameter \code{top}.
#'
#' \code{model.matrix(which = "random", top = NULL)} returns the result
#' of combining the parameter-specific all-terms matrices \emph{and}
#' permuting the columns to obtain a convenient ordering of random effects
#' coefficients.
#' (Coefficients are sorted by relation to a common random vector;
#' random vectors are sorted by relation to a common covariance matrix.)
#'
#' None of these combined matrices are "proper" design matrices,
#' as none possess an \code{assign} attribute.
#'
#' @return
#' A \link[=matrix]{dense} or \link[Matrix:sparseMatrix]{sparse} matrix
#' with attributes \code{assign} and \code{contrasts}.
#' \code{assign} is absent in special cases; see Details.
#'
#' @family extractors
#' @method model.matrix egf
#' @export
#' @importFrom stats formula terms model.frame model.matrix
model.matrix.egf <- function(object,
                             which = c("fixed", "random"),
                             top = NULL,
                             random = NULL,
                             ...) {
  which <- match.arg(which)

  if (is.null(top)) {
    ## Return the combined fixed or random effects design matrix
    ## from object internals
    name <- switch(which, fixed = if (object$control$sparse_X) "Xs" else "Xd", random = "Z")
    res <- object$tmb_out$env$data[[name]]

    ## Append 'contrasts' but not 'assign', which only makes sense
    ## for submatrices
    attr(res, "contrasts") <- object$contrasts[[substr(name, 1L, 1L)]]
    return(res)
  }

  top <- match.arg(top, egf_get_names_top(object, link = TRUE))
  frame <- model.frame(object, which = "parameters", top = top)
  l <- split_effects(formula(terms(frame))) # list(fixed = <formula>, random = <list of calls to `|`>)

  if (which == "fixed") {
    ## Return parameter-specific fixed effects design matrix
    res <- egf_make_X(fixed = l$fixed, data = frame, sparse = object$control$sparse_X)
    return(res)
  }

  random <- substitute(random)
  any_random_effects <- length(l$random) > 0L

  if (is.null(random)) {
    if (!any_random_effects) {
      ## Return empty sparse matrix with correct number of rows
      res <- object$tmb_out$env$data$Z[, integer(0L), drop = FALSE]
      return(res)
    }

    ## Return parameter-specific combined random effects design matrix
    Z <- lapply(l$random, egf_make_Z, data = frame)
    res <- do.call(cbind, Z)

    ## Append 'contrasts' but not 'assign', which only makes sense
    ## for submatrices
    contrasts <- unlist1(lapply(Z, attr, "contrasts"))
    attr(res, "contrasts") <- contrasts[!duplicated(names(contrasts))]
    return(res)
  }

  if (!any_random_effects) {
    stop(wrap("Expected 'random = NULL': mixed effects formula for parameter ", sQuote(top), " does not contain random effects terms."))
  }

  l$random[] <- lapply(l$random, function(x) call("(", x))
  if (!any(l$random == random)) {
    stop(wrap("Expected 'random = NULL' or 'random' matching one of:"), "\n\n", paste0("  ", l$random, collapse = "\n"))
  }

  ## Return term-specific random effects design matrix
  egf_make_Z(random = random[[2L]], data = frame)
}

#' @rdname model.matrix.egf
#' @method model.matrix egf_no_fit
#' @export
#' @importFrom stats formula terms model.frame model.matrix
model.matrix.egf_no_fit <- model.matrix.egf

#' Extract model terms
#'
#' Extract \code{"\link[=terms.object]{terms}"} objects
#' corresponding to top level nonlinear model parameters.
#'
#' @param x
#'   An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.
#' @param top
#'   A character string specifying a top level nonlinear model parameter.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A \code{"\link[=terms.object]{terms}"} object.
#'
#' @family extractors
#' @export
#' @importFrom stats terms model.frame
terms.egf <- function(x, top = egf_get_names_top(x, link = TRUE), ...) {
  top <- match.arg(top)
  frame <- model.frame(x, which = "parameters", top = top)
  terms(frame)
}

#' @rdname terms.egf
#' @export
#' @importFrom stats terms model.frame
terms.egf_no_fit <- terms.egf

#' Extract model formulae
#'
#' Extract mixed effects model formulae corresponding to top level
#' nonlinear model parameters.
#'
#' @param x
#'   An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.
#' @param top
#'   A character string specifying a top level nonlinear model parameter.
#' @param split
#'   A logical flag. If \code{TRUE}, then fixed and random effects terms
#'   are returned separately.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' If \code{split = FALSE}, then the mixed effects model formula
#' corresponding to parameter \code{top}.
#' If \code{split = TRUE}, then the same formula but with random effects
#' terms deleted from the right hand side.
#' The deleted terms are preserved in a list and retained as an attribute
#' of the result, namely \code{random}.
#'
#' @family extractors
#' @export
#' @importFrom stats formula terms
formula.egf <- function(x, top = egf_get_names_top(x, link = TRUE), split = FALSE, ...) {
  top <- match.arg(top)
  res <- formula(terms(x, top = top))
  if (split) {
    l <- split_effects(res)
    res <- l$fixed
    attr(res, "random") <- lapply(l$random, function(x) call("(", x))
  }
  res
}

#' @rdname formula.egf
#' @export
#' @importFrom stats formula terms
formula.egf_no_fit <- formula.egf

#' Extract number of observations
#'
#' Returns the number of observations of disease incidence that were used
#' or would be used in estimation of a model. This number excludes missing
#' values and observations not belonging to a fitting window, which,
#' despite being preserved in model objects, do not affect estimation.
#'
#' @param object
#'   An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' An integer.
#'
#' @family extractors
#' @export
#' @importFrom stats model.frame nobs
nobs.egf <- function(object, ...) {
  sum(!is.na(model.frame(object)$x))
}

#' @rdname nobs.egf
#' @export
#' @importFrom stats model.frame nobs
nobs.egf_no_fit <- nobs.egf

#' Extract residual degrees of freedom
#'
#' Returns the number of observations (see \code{\link[=nobs.egf]{nobs}})
#' minus the number of estimated parameters (fixed effects coefficients
#' and random effect covariance parameters).
#'
#' @param object
#'   An \code{"\link{egf}"} or \code{"\link[=egf]{egf_no_fit}"} object.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' An integer.
#'
#' @family extractors
#' @method df.residual egf
#' @export
#' @importFrom stats nobs df.residual
df.residual.egf <- function(object, ...) {
  as.numeric(nobs(object)) - sum(!object$random)
}

#' @rdname df.residual.egf
#' @export
#' @importFrom stats nobs df.residual
df.residual.egf_no_fit <- df.residual.egf

#' Extract log likelihood
#'
#' Retrieves the value of log likelihood from a fitted model object.
#'
#' @param object
#'   An \code{"\link{egf}"} object.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A number inheriting from class \code{"logLik"}.
#' Attribute \code{df} is the number of estimated parameters
#' (fixed effects coefficients and random effect covariance parameters).
#' Attribute \code{nobs} is the number of observations of disease incidence
#' used in estimation.
#'
#' @family extractors
#' @export
#' @importFrom stats nobs logLik
logLik.egf <- function(object, ...) {
  res <- -object$value
  attr(res, "df") <- sum(!object$random)
  attr(res, "nobs") <- nobs(object)
  class(res) <- "logLik"
  res
}

#' @export
#' @importFrom stats extractAIC
extractAIC.egf <- function(fit, scale, k = 2, ...) {
  ll <- logLik(fit)
  edf <- attr(ll, "df")
  c(edf, -2 * as.numeric(ll) + k * edf)
}
