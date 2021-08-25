#' Construct data objects for C++ template
#'
#' Gathers in a \link{list} data objects to be passed to
#' the package's C++ template via \pkg{TMB}'s \code{DATA_*} macros.
#'
#' @inheritParams egf
#' @param frame,frame_parameters
#'   Model frames obtained from the list output of \code{egf_make_frames}.
#'
#' @return
#' [Below,
#' \code{N = \link{nlevels}(frame$window)}
#' is the number of fitting windows,
#' \code{n = N + \link{sum}(!\link{is.na}(frame$window))}
#' is the total number of time points associated with a fitting window, and
#' \code{p = \link{length}(frame_parameters)}
#' is the number of top level nonlinear model parameters.]
#'
#' A \link{list} inheriting from \link{class} \code{"tmb_data"},
#' with elements:
#' \item{time}{
#'   A \link[=double]{numeric} vector of length \code{n} giving
#'   times since the left endpoint of the current fitting window.
#' }
#' \item{time_seg_len}{
#'   An \link{integer} vector of length \code{N} specifying the
#'   length of each fitting window as a number of time points.
#' }
#' \item{x}{
#'   A \link[=double]{numeric} vector of length \code{n-N} giving incidence
#'   in each fitting window. \code{x[i]} in window \code{k} is the number
#'   of cases observed from \code{time[k+i-1]} to \code{time[k+i]}.
#' }
#' \item{day1}{
#'   If \code{model$day_of_week > 0}, then an \link{integer} vector
#'   of length \code{N} indicating the first day of week in each
#'   fitting window, with value \code{i} in \code{0:6} mapping to
#'   the day of week \code{i} days after the reference day specified
#'   by \code{model$day_of_week}. Otherwise, an integer vector of
#'   the same length filled with \code{-1}.
#' }
#' \item{flags}{
#'   A \link{list} with \link{integer} elements, used as flags to
#'   specify the model being estimated and to indicate what blocks
#'   of template code should be run.
#' }
#' \item{indices}{
#'   A \link{list} with \link{integer} elements and \link{names} of the
#'   form \code{"index_link_parameter"} (e.g., \code{"index_log_r"}),
#'   giving the column 0-index of top level nonlinear model parameters
#'   (e.g., \code{log(r)}) in the response matrix. Value \code{-1}
#'   is used for parameters not belonging to the model being estimated.
#' }
#' \item{Y}{
#'   The \link[=model.offset]{offset} component of the response matrix
#'   in dense format, with \code{N} rows and \code{p} columns.
#' }
#' \item{X}{
#'   The fixed effects design matrix in \link[Matrix:sparseMatrix]{sparse}
#'   or \link[=matrix]{dense} format (depending on \code{control$sparse_X}),
#'   with \code{N} rows.
#' }
#' \item{Xs, Xd}{
#'   If \code{control$sparse_X = TRUE}, then \code{Xs = X} and \code{Xd}
#'   is an empty dense matrix. Otherwise, \code{Xd = X} and \code{Xs} is
#'   an empty sparse matrix.
#' }
#' \item{Z}{
#'   The random effects design matrix in \link[Matrix:sparseMatrix]{sparse}
#'   format, with \code{N} rows. If there are no random effects, then \code{Z}
#'   is an empty sparse matrix.
#' }
#' \item{beta_index, b_index}{
#'   \link[=integer]{Integer} vectors of length \code{\link{ncol}(X)}
#'   and \code{\link{ncol}(Z)}, respectively, with values in \code{0:(p-1)}.
#'   These split the columns of \code{X} and \code{Z} by relation to
#'   a common top level nonlinear model parameter.
#' }
#' \item{beta_index_tab, b_index_tab}{
#'   \link[=integer]{Integer} vectors of length \code{p} counting the
#'   columns of \code{X} and \code{Z}, respectively, that relate to
#'   a common top level nonlinear model parameter.
#' }
#' \item{block_rows, block_cols}{
#'   \link[=integer]{Integer} vectors together giving the dimensions
#'   of each block of random effects.
#' }
#'
#' @seealso \code{\link[TMB]{MakeADFun}}
#' @noRd
#' @importFrom stats formula terms model.offset
#' @importFrom Matrix sparseMatrix
egf_make_tmb_data <- function(model, frame, frame_parameters, control, fit, init) {
  ## Indices of time points associated with fitting windows
  first <- attr(frame, "first")
  last <- attr(frame, "last")
  index <- Map(seq.int, first, last)

  ## Fitting window lengths as numbers of time points
  time_seg_len <- lengths(index)
  N <- length(index)

  ## Time since earliest time point
  ulindex <- unlist(index, FALSE, FALSE)
  time <- frame$time[ulindex] - rep.int(frame$time[first], time_seg_len)

  ## Incidence
  x <- frame$x[!is.na(frame$window)]

  if (model$day_of_week > 0L) {
    ## Date during 24 hours starting at earliest time point
    Date1 <- .Date(frame$time[first])

    ## Day of week on that Date coded as an integer `i` in `0:6`.
    ## Integer `i` maps to the day of week `i` days after a reference day,
    ## which is the day `model$day_of_week` days after Saturday
    ## NB: weekdays(.Date(2L)) == "Saturday"
    day1 <- as.integer(julian(Date1, origin = .Date(2L + model$day_of_week)) %% 7)
  } else {
    day1 <- rep_len(-1L, N)
  }

  ## Response matrix, offset component only
  offsets <- lapply(frame_parameters, model.offset)
  offsets[vapply(offsets, is.null, FALSE)] <- list(rep_len(0, N))
  Y <- do.call(cbind, offsets)

  ## Top level nonlinear model parameter names
  names_top <- names(frame_parameters)

  ## List of fixed effects formulae and list of lists of random effects terms
  formula_parameters <- lapply(frame_parameters, function(x) formula(terms(x)))
  effects <- lapply(formula_parameters, split_effects)
  fixed <- lapply(effects, `[[`, "fixed")
  random <- lapply(effects, `[[`, "random")

  ## Fixed effects infrastructure
  X <- Map(egf_make_X, fixed = fixed, data = frame_parameters, sparse = control$sparse_X)
  X <- egf_combine_X(fixed = fixed, X = X)
  X_info <- attr(X, "info")
  beta_index <- as.integer(X_info$top) - 1L
  beta_index_tab <- c(table(X_info$top))

  ## Random effects infrastructure
  random1 <- do.call(c, unname(random))
  names(random1) <- rep.int(names(random), lengths(random))
  Z <- Map(egf_make_Z, random = random1, data = rep.int(frame_parameters, lengths(random)))
  Z <- egf_combine_Z(random = random1, Z = Z)
  if (is.null(Z)) {
    b_index <- integer(0L)
    b_index_tab <- rep_len(0L, length(names_top))
    block_rows <- integer(0L)
    block_cols <- integer(0L)
  } else {
    attr(Z, "info")$top <- factor(attr(Z, "info")$top, levels = names_top)
    Z_info <- attr(Z, "info")
    b_index <- as.integer(Z_info$top) - 1L
    b_index_tab <- c(table(Z_info$top))
    cor <- interaction(Z_info[c("term", "group")], drop = TRUE, lex.order = TRUE)
    vec <- interaction(Z_info[c("term", "group", "level")], drop = TRUE, lex.order = TRUE)
    block_rows <- as.integer(colSums(table(Z_info$top, cor) > 0L))
    block_cols <- as.integer(colSums(table(vec, cor) > 0L))
  }

  ## Flags
  flags <- list(
    flag_curve = get_flag("curve", model$curve),
    flag_excess = as.integer(model$excess),
    flag_family = get_flag("family", model$family),
    flag_day_of_week = as.integer(model$day_of_week > 0L),
    flag_trace = control$trace,
    flag_sparse_X = as.integer(control$sparse_X),
    flag_predict = 0L
  )

  ## Column indices of top level nonlinear model parameters
  ## in response matrix
  names_top_all <- egf_get_names_top(NULL, link = TRUE)
  indices <- as.list(match(names_top_all, names_top, 0L) - 1L)
  names(indices) <- sub("^(log|logit)\\((.*)\\)$", "index_\\1_\\2", names_top_all)

  ## Empty design matrices
  empty_dense_matrix <- `dim<-`(integer(0L), c(N, 0L))
  empty_sparse_matrix <- sparseMatrix(
    i = integer(0L),
    j = integer(0L),
    x = integer(0L),
    dims = c(N, 0L)
  )

  res <- list(
    time = time,
    time_seg_len = time_seg_len,
    x = x,
    day1 = day1,
    flags = flags,
    indices = indices,
    Y = Y,
    X = X,
    Xd = if (control$sparse_X) empty_dense_matrix else X,
    Xs = if (control$sparse_X) X else empty_sparse_matrix,
    Z = if (is.null(Z)) empty_sparse_matrix else Z,
    beta_index = beta_index,
    b_index = b_index,
    beta_index_tab = beta_index_tab,
    b_index_tab = b_index_tab,
    block_rows = block_rows,
    block_cols = block_cols
  )
  class(res) <- c("tmb_data", "list")
  res
}

#' Construct parameter objects for C++ template
#'
#' Gathers in a \link{list} parameter objects \code{beta}, \code{theta},
#' and \code{b} to be passed to the package's C++ template via \pkg{TMB}'s
#' \code{PARAMETER_*} macros during the first likelihood evaluation.
#' See also \code{\link[TMB]{MakeADFun}}.
#'
#' @inheritParams egf
#' @param tmb_data
#'   A \code{"\link[=egf_make_tmb_data]{tmb_data}"} object.
#' @param
#'   Model frames obtained from the list output of \code{egf_make_frames}.
#'
#' @details
#' When \code{init = NULL}, naive estimates of top level nonlinear
#' model parameters are obtained for each fitting window as follows:
#' \describe{
#' \item{\code{r}}{
#'   The slope of a linear model fit to \code{\link{log1p}(\link{cumsum}(x)))}.
#' }
#' \item{\code{alpha}}{
#'   \code{r*c0^(1-p)} if \code{curve = "subexponential"},
#'   \code{r/log(K/c0)} if \code{curve = "gompertz"}.
#'   These are the values obtained by setting the per capita growth
#'   rate at time 0 in the subexponential and Gompertz models equal
#'   to \code{r}, substituting the naive estimates of \code{r},
#'   \code{c0}, \code{K}, and \code{p}, and solving for \code{alpha}.
#' }
#' \item{\code{c0}}{
#'   \code{\link{exp}(log_c0)}, where \code{log_c0} is the intercept
#'   of a linear model fit to \code{\link{log1p}(\link{cumsum}(x))}.
#' }
#' \item{\code{tinfl}}{
#'   \code{\link{max}(time)}. This assumes that the fitting window
#'   ends near the time of a peak in incidence.
#' }
#' \item{\code{K}}{
#'   \code{2*\link{sum}(x)}. This assumes that the fitting window
#'   ends near the time of a peak in incidence _and_ that incidence
#'   is roughly symmetric about the peak.
#' }
#' \item{\code{p}}{0.95}
#' \item{\code{a, b, disp, w[123456]}}{1}
#' }
#' The naive estimates are log- or logit-transformed (all but \code{p}
#' use a log link), and the initial value of the \code{"(Intercept)"}
#' coefficient in a top level parameter's fixed effects model
#' (if that coefficient exists) is taken to be the mean of the link scale
#' estimates across fitting windows. The remaining elements of \code{beta}
#' take initial value 0. All elements of \code{theta} and \code{b}
#' take initial value 0, so that, initially, all random vectors follow
#' a standard normal distribution.
#'
#' @return
#' A \link{list} with elements:
#' \item{beta}{
#'   A \link[=double]{numeric} vector of length
#'   \code{\link{length}(tmb_data$beta_index)} listing
#'   initial values for fixed effects coefficients.
#' }
#' \item{theta}{
#'   A \link[=double]{numeric} vector of length
#'   \code{\link{sum}(f(tmb_data$block_rows))}, where \code{f(n) = n*(n+1)/2},
#'   specifying initial covariance matrices for all random vectors.
#'   Let \code{x} be a segment of \code{theta} of length \code{f(n)},
#'   corresponding to a random vector \code{u} of length \code{n}.
#'   Then \code{x[1:n]} lists log standard deviations of the elements
#'   of \code{u}, and \code{x[-(1:n)]} lists (in row-major order)
#'   the subdiagonal elements of the unit lower triangular matrix
#'   \code{L} in the Cholesky factorization of the correlation matrix
#'   \code{S} of \code{u}:\cr
#'   \code{S = A \link[=matmult]{\%*\%} \link{t}(A)}\cr
#'   \code{A = D^-0.5 \link[=matmult]{\%*\%} L}\cr
#'   \code{D = \link{diag}(diag(L \link[=matmult]{\%*\%} \link{t}(L)))}
#' }
#' \item{b}{
#'   A \link[=double]{numeric} vector of length
#'   \code{\link{length}(tmb_data$b_index)} listing initial values
#'   for random effects coefficients (unit variance scale).
#' }
#'
#' @noRd
#' @importFrom stats coef lm na.omit qlogis terms
egf_make_tmb_parameters <- function(tmb_data, model, frame, frame_parameters,
                                    fit, init) {
  ## Lengths of parameter objects
  f <- function(n) as.integer(n * (n + 1) / 2)
  len <- c(
    beta  = length(tmb_data$beta_index),
    theta = sum(f(tmb_data$block_rows)),
    b     = length(tmb_data$b_index)
  )

  ## If user does not specify a full parameter vector
  if (is.null(init)) {
    ## Initialize each parameter object to a vector of zeros
    init_split <- lapply(len, numeric)

    ## Identify top level nonlinear model parameters
    ## whose mixed effects formulae have an intercept
    has_intercept <- function(frame) {
      attr(terms(frame), "intercept") == 1L
    }
    i1 <- vapply(frame_parameters, has_intercept, FALSE)

    ## For each of these top level nonlinear model parameters,
    ## compute the mean over all fitting windows of the naive estimate,
    ## and assign the result to the coefficient of `beta` corresponding
    ## to "(Intercept)"
    if (any(i1)) {
      ## Time series segments
      tx_split <- split(frame[c("time", "x")], frame$window)

      ## Functions computing naive estimates given time series segments
      get_r_c0 <- function(d) {
        n <- max(2, trunc(nrow(d) / 2))
        ab <- try(coef(lm(log1p(cumsum(x)) ~ time, data = d, subset = seq_len(n), na.action = na.omit)), silent = TRUE)
        if (inherits(ab, "try-error") || !all(is.finite(ab))) {
          return(c(0.04, 1))
        }
        c(ab[[2L]], exp(ab[[1L]]))
      }
      get_tinfl <- function(d) {
        max(d$t)
      }
      get_K <- function(d) {
        2 * sum(d$x, na.rm = TRUE)
      }

      ## Naive estimates for each fitting window
      r_c0 <- vapply(tx_split, get_r_c0, c(0, 0))
      r  <- r_c0[1L, ]
      c0 <- r_c0[2L, ]
      tinfl <- vapply(tx_split, get_tinfl, 0)
      K <- vapply(tx_split, get_K, 0)
      p <- 0.95
      alpha <- switch(model$curve,
        subexponential = r * c0^(1 - p),
        gompertz = r / (log(K) - log(c0)),
        NA_real_
      )
      Y_init <- data.frame(r, alpha, c0, tinfl, K, p)
      Y_init[c("a", "b", "disp", paste0("w", 1:6))] <- 1

      ## Link transform
      Y_init[] <- Map(function(x, s) egf_link_match(s)(x),
        x = Y_init,
        s = egf_link_get(names(Y_init))
      )
      names(Y_init) <- egf_link_add(names(Y_init))

      ## Identify elements of `beta` corresponding to "(Intercept)"
      ## and assign means over fitting windows
      nfp <- names(frame_parameters)[i1]
      index <- match(nfp, attr(tmb_data$X, "info")$top, 0L)
      init_split$beta[index] <- colMeans(Y_init[nfp])
    }
  } else {
    ## Validate and split full parameter vector,
    ## producing `list(beta, b, theta)`
    stop_if_not(
      is.numeric(init),
      length(init) == sum(len),
      is.finite(init),
      m = sprintf("`init` must be a finite numeric vector of length %d.", sum(len))
    )
    names(init) <- NULL
    init_split <- split(init, rep.int(gl(3L, 1L, labels = names(len)), len))
  }

  ## Retain all naive estimates if debugging
  if (is.null(init) && !fit) {
    attr(init_split, "Y_init") <- as.matrix(Y_init[names(frame_parameters)])
  }
  init_split
}

egf_make_tmb_args <- function(model, frame, frame_parameters, control, fit, init, map) {
  tmb_data <- egf_make_tmb_data(
    model = model,
    frame = frame,
    frame_parameters = frame_parameters,
    control = control
  )
  tmb_parameters <- egf_make_tmb_parameters(
    tmb_data = tmb_data,
    model = model,
    frame = frame,
    frame_parameters = frame_parameters,
    fit = fit,
    init = init
  )
  if (is.null(map)) {
    tmb_map <- list()
  } else {
    len <- lengths(tmb_parameters)
    stop_if_not(
      is.factor(map),
      length(map) == sum(len),
      m = sprintf("`map` must be a factor of length %d.", sum(len))
    )
    tmb_map <- split(map, rep.int(gl(3L, 1L, labels = names(len)), len))
    tmb_map <- lapply(tmb_map, factor)
  }
  if (egf_has_random(tmb_data)) {
    ## Declare that `b` contains random effects
    tmb_random <- "b"
  } else {
    ## Declare that there are no random effects
    tmb_random <- NULL
    ## Fix `theta` and `b` to NA_real_ since only `beta` is used
    tmb_parameters$theta <- tmb_parameters$b <- NA_real_
    tmb_map$theta <- tmb_map$b <- factor(NA)
  }
  list(
    data = tmb_data,
    parameters = tmb_parameters,
    map = tmb_map,
    random = tmb_random,
    profile = if (control$profile) "beta" else NULL,
    DLL = "epigrowthfit",
    silent = (control$trace == 0L)
  )
}

egf_update_tmb_args <- function(tmb_args, priors_top, priors_bottom) {
  f <- function(priors) {
    n <- length(priors)
    has_prior <- !vapply(priors, is.null, FALSE)
    flag_regularize <- rep_len(-1L, n)
    flag_regularize[has_prior] <- get_flag("prior", vapply(priors[has_prior], `[[`, "", "family"))
    hyperparameters <- rep_len(list(numeric(0L)), n)
    hyperparameters[has_prior] <- lapply(priors[has_prior], function(x) unlist(x$parameters, FALSE, FALSE))
    list(flag_regularize = flag_regularize, hyperparameters = hyperparameters)
  }
  l1 <- f(priors_top)
  tmb_args$data$flags$flag_regularize_top <- l1$flag_regularize
  tmb_args$data$hyperparameters_top <- l1$hyperparameters
  l2 <- f(priors_bottom)
  tmb_args$data$flags$flag_regularize_bottom <- l2$flag_regularize
  tmb_args$data$hyperparameters_bottom <- l2$hyperparameters
  tmb_args
}

#' Patch TMB-generated functions
#'
#' Define wrapper functions on top of \code{\link[TMB]{MakeADFun}}-generated
#' functions \code{fn} and \code{gr}, so that function and gradient evaluations
#' can retry inner optimization using fallback methods in the event that the
#' default method (usually \code{\link[TMB]{newton}}) fails.
#'
#' @param fn,gr
#'   \link[=function]{Function}s from a \code{\link[TMB]{MakeADFun}}-generated
#'   \link{list} object to be patched.
#' @param inner_optimizer
#'   A \link{list} of \code{"\link{egf_inner_optimizer}"} objects
#'   specifying inner optimization methods to be tried in turn.
#'
#' @return
#' A patched version of \link{function} \code{fn} or \code{gr}.
#'
#' @noRd
NULL

egf_patch_fn <- function(fn, inner_optimizer) {
  e <- environment(fn)
  if (!exists(".egf_env", where = e, mode = "environment", inherits = FALSE)) {
    e$.egf_env <- new.env(parent = emptyenv())
  }
  e$.egf_env$fn <- fn
  e$.egf_env$inner_optimizer <- inner_optimizer

  last.par <- random <- inner.method <- inner.control <- .egf_env <- NULL # for `check`
  pfn <- function(x = last.par[-random], ...) {
    oim <- inner.method
    oic <- inner.control
    on.exit({
      inner.method <<- oim
      inner.control <<- oic
    })
    for (io in .egf_env$inner_optimizer) {
      inner.method <<- io$method
      inner.control <<- io$control
      v <- .egf_env$fn(x, ...)
      if (is.numeric(v) && length(v) == 1L && is.finite(v)) {
        return(v)
      }
    }
    NaN # no warning to avoid duplication of `optim` and `nlminb` warnings
  }
  environment(pfn) <- e
  pfn
}

egf_patch_gr <- function(gr, inner_optimizer) {
  e <- environment(gr)
  if (!exists(".egf_env", where = e, mode = "environment", inherits = FALSE)) {
    e$.egf_env <- new.env(parent = emptyenv())
  }
  e$.egf_env$gr <- gr
  e$.egf_env$inner_optimizer <- inner_optimizer

  last.par <- random <- inner.method <- inner.control <- .egf_env <- NULL # for `check`
  pgr <- function(x = last.par[-random], ...) {
    oim <- inner.method
    oic <- inner.control
    on.exit({
      inner.method <<- oim
      inner.control <<- oic
    })
    n <- length(x)
    for (io in .egf_env$inner_optimizer) {
      inner.method <<- io$method
      inner.control <<- io$control
      v <- .egf_env$gr(x, ...)
      if (is.numeric(v) && length(v) == n && all(is.finite(v))) {
        return(v)
      }
    }
    warning("Unable to evaluate `gr(x)`, returning NaN.")
    NaN # warning because scalar result is unexpected
  }
  environment(pgr) <- e
  pgr
}
