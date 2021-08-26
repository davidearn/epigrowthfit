#' Simulate incidence from a fitted model
#'
#' Simulates new incidence data conditional on a fitted nonlinear
#' mixed effects model of epidemic growth. Only observations within
#' fitting windows are simulated.
#'
#' @param object
#'   An \code{"\link{egf}"} object specifying a fitted nonlinear
#'   mixed effects model.
#' @param nsim
#'   A positive integer indicating a number of replications.
#' @param seed
#'   An integer used to set the \link{RNG} state before simulation.
#'   The default (\code{\link{NULL}}) is to use the state at the
#'   time of the function call. The RNG state
#'   (either a \link{list} of arguments to \code{\link{set.seed}}
#'   or a value of \code{\link{.Random.seed}}) is preserved as an
#'   \link[=attributes]{attribute} of the result.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A \link[=data.frame]{data frame} inheriting from \link{class}
#' \code{"egf_simulate"}, with variables \code{ts}, \code{window},
#' and \code{time}, and \code{nsim} further variables with names
#' of the form \code{x[0-9]+}.
#' It corresponds rowwise to \code{object$frame[!is.na(object$frame$window), ]}.
#'
#' @export
#' @importFrom stats simulate
simulate.egf <- function(object, nsim = 1L, seed = NULL, ...) {
  stop_if_not_integer(nsim, "positive")
  nsim <- as.integer(nsim)
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    runif(1L)
  }
  if (is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    oRNGstate <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", oRNGstate, envir = .GlobalEnv))
    RNGstate <- c(list(seed), as.list(RNGkind()))
    names(RNGstate) <- names(formals(set.seed))
    do.call(set.seed, RNGstate)
  }

  frame <- object$frame[!is.na(object$frame$window), c("ts", "window", "time"), drop = FALSE]
  nx <- sprintf("x%0*d", as.integer(log10(nsim)) + 1L, seq_len(nsim))
  frame[nx] <- replicate(nsim, object$tmb_out$simulate()$x)
  attr(frame, "RNGstate") <- RNGstate
  class(frame) <- c("egf_simulate", "data.frame")
  frame
}

#' Simulate incidence time series
#'
#' Simulates incidence time series with daily observations according to
#' a specified nonlinear model. Top level nonlinear model parameters vary
#' between time series according to a fixed intercept model \code{~ts} or
#' random intercept model \code{~(1 | ts)}.
#'
#' @param object
#'   An \code{"\link{egf_model}"} specifying a top level nonlinear model
#'   to be simulated.
#' @param nsim
#'   A positive integer indicating a number of simulated time series.
#' @param seed
#'   An integer used to set the \link{RNG} state before simulation.
#'   The default (\code{\link{NULL}}) is to use the state at the
#'   time of the function call. The RNG state
#'   (either a \link{list} of arguments to \code{\link{set.seed}}
#'   or a value of \code{\link{.Random.seed}}) is preserved as an
#'   \link[=attributes]{attribute} of the result.
#' @param mu
#'   A \link{numeric} vector listing means across time series
#'   of top level nonlinear model parameters (link scale).
#'   It is assumed that elements are ordered as in
#'   \code{\link{egf_get_names_top}(object, link = TRUE)}.
#' @param Sigma
#'   A symmetric positive definite \link{numeric} \link{matrix}
#'   to be used as the covariance matrix corresponding to \code{mu}.
#'   The default (\code{\link{NULL}}) is equivalent (conceptually)
#'   to a zero matrix and is handled specially.
#' @param tol
#'   A non-negative number indicating a tolerance for lack of positive
#'   definiteness of \code{Sigma}. Negative eigenvalues of \code{Sigma}
#'   must not be less than \code{-tol * rho}, where \code{rho} is the
#'   spectral radius of \code{Sigma}. (However, regardless of \code{tol},
#'   \code{\link{diag}(Sigma)} must be positive, as standard deviations
#'   are processed on the log scale.)
#' @param tmax
#'   A positive number. Simulated time series run from 0 days to 1 day
#'   after the inflection time of the expected cumulative incidence
#'   function indicated by \code{object$curve}. However, if there is no
#'   inflection time, as in the exponential and subexponential models
#'   of cumulative incidence, then simulated time series run from 0 days
#'   to \code{tmax} days.
#' @param cstart
#'   A number indicating a threshold value of cumulative incidence.
#'   Left endpoints of suggested fitting windows are those times
#'   when cumulative incidence first exceeds this threshold.
#' @param ...
#'   Unused optional arguments.
#'
#' @return
#' A \link{list} inheriting from \link{class} \code{"egf_model_simulate"},
#' with elements:
#' \item{model, mu, Sigma, origin}{
#'   Copies of the so-named arguments.
#'   (\code{model} is a copy of argument \code{object}.)
#' }
#' \item{Y}{
#'   A \link[=double]{numeric} \link{matrix} with \code{nsim} rows and
#'   \code{\link{length}(mu)} columns listing the top level nonlinear model
#'   parameter values underlying each time series.
#'   If \code{Sigma} is \code{\link{NULL}}, then the row vectors of \code{Y}
#'   are all \code{mu}. Otherwise, \code{Y} is (conceptually) the result of
#'   \code{MASS::mvrnorm(nsim, mu, Sigma, tol)}.
#' }
#' \item{formula}{
#'   A \link{formula}, \code{cbind(time, x) ~ ts}, expressing how
#'   simulated time series are stored in \code{data}.
#' }
#' \item{formula_windows}{
#'   A \link{formula}, \code{cbind(start, end) ~ ts}, expressing how
#'   fitting window endpoints are stored in \code{data_windows}.
#' }
#' \item{formula_parameters}{
#'   A \link{formula} specifying the generative model.
#'   If \code{Sigma = \link{NULL}}, then the formula is
#'   \code{~1} if \code{nsim = 1} and \code{~ts} if \code{nsim > 1}.
#'   Otherwise, it is \code{~(1 | ts)}.
#' }
#' \item{data}{
#'   A \link[=data.frame]{data frame} with variables \code{ts}, \code{time},
#'   and \code{x} storing \code{nsim} simulated time series in long format.
#' }
#' \item{data_windows}{
#'   A \link[=data.frame]{data frame} with \code{nsim} rows and variables
#'   \code{ts}, \code{start}, and \code{end} suggesting a fitting window
#'   for each simulated time series. Start times are determined according
#'   to \code{cstart}. End times are precisely the last time point in the
#'   corresponding time series.
#' }
#' \item{actual}{
#'   A \link[=double]{numeric} vector giving the full parameter vector
#'   of the generative model. When estimating this model from \code{data},
#'   \code{\link{egf}} output should be compared against \code{actual}.
#'   More precisely, if \code{m} is the \code{"\link{egf}"} object,
#'   then \code{m$best} estimates \code{actual}.
#' }
#' \item{nonrandom}{
#'   An \link{integer} vector indexing the elements of \code{actual}
#'   that are not random effects.
#' }
#' \item{call}{
#'   The \link{call} to \code{simulate.egf_model}, allowing for updates
#'   to the \code{"egf_model_simulate"} object via \code{\link{update}}.
#' }
#'
#' @examples
#' model <- egf_model(curve = "logistic", family = "nbinom")
#'
#' r <- log(2) / 20
#' tinfl <- 160
#' K <- 25000
#' disp <- 50
#'
#' mu <- log(c(r, tinfl, K, disp))
#' Sigma <- diag(rep_len(0.5^2, length(mu)))
#'
#' sim <- simulate(model,
#'   nsim = 20L,
#'   seed = 202737L,
#'   mu = mu,
#'   Sigma = Sigma,
#'   cstart = 10
#' )
#'
#' @seealso
#' \code{\link{egf.egf_model_simulate}} for estimating the generative model
#' from simulated time series,
#' \code{\link{simulate.egf}} for simulating time series from \emph{fitted}
#' models
#'
#' @export
#' @importFrom stats simulate update runif cov2cor
simulate.egf_model <- function(object, nsim = 1L, seed = NULL,
                               mu, Sigma = NULL, tol = 1e-06,
                               tmax = 100, cstart = 0, ...) {
  stop_if_not_integer(nsim, "positive")
  nsim <- as.integer(nsim)
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    runif(1L)
  }
  if (is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    set_RNGstate <- function() assign(".Random.seed", RNGstate, .GlobalEnv)
  } else {
    oRNGstate <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", oRNGstate, envir = .GlobalEnv))
    RNGstate <- c(list(as.numeric(seed)), as.list(RNGkind()))
    names(RNGstate) <- names(formals(set.seed))
    set_RNGstate <- function() do.call(set.seed, RNGstate)
  }

  names_top <- egf_get_names_top(object, link = TRUE)
  p <- length(names_top)
  stopifnot(
    is.numeric(mu),
    length(mu) == p,
    is.finite(mu)
  )
  names(mu) <- names_top
  if (!is.null(Sigma)) {
    stopifnot(
      is.matrix(Sigma),
      is.numeric(Sigma),
      dim(Sigma) == length(mu),
      is.finite(Sigma),
      isSymmetric(Sigma)
    )
    stop_if_not_number(tol, "nonnegative")
    e <- eigen(Sigma, symmetric = TRUE)$values
    stop_if_not(
      e >= -tol * abs(e[1L]),
      m = "'Sigma' must be positive definite."
    )
    dimnames(Sigma) <- rep_len(list(names_top), 2L)
  }

  has_inflection <- object$curve %in% c("gompertz", "logistic", "richards")
  if (has_inflection) {
    ## Time: daily from 0 days to (inflection time)+1 days
    if (is.null(Sigma)) {
      ## Fixed intercept model: inflection time is known up front.
      tmax <- ceiling(exp(mu[["log(tinfl)"]])) + 1
    } else {
      ## Random intercept model: inflection time is not known up front.
      ## Two passes are necessary. The first pass asks for time series
      ## of minimal length (hence 'tmax <- 1') and serves only to
      ## retrieve randomly generated inflection times. The second pass
      ## asks for time series of appropriate length.
      tmax <- 1
    }
  } else {
    ## Time: daily from 0 days to user-specified 'tmax' days
    stop_if_not_number(tmax, "positive")
    tmax <- max(1, trunc(tmax))
  }
  stop_if_not_number(cstart)

  ## Construct arguments to 'egf' corresponding to 'nsim', 'mu', 'Sigma'
  formula <- cbind(time, x) ~ ts
  formula_windows <- cbind(start, end) ~ ts
  data <- data.frame(
    ts = gl(nsim, 1L + as.integer(tmax)),
    time = seq.int(0, tmax, by = 1),
    x = 0L # arbitrary non-negative integer to pass checks
  )
  data_windows <- data.frame(
    ts = gl(nsim, 1L),
    start = 0,
    end = tmax
  )
  if (is.null(Sigma)) {
    if (nsim == 1L) {
      formula_parameters <- ~1
      init <- mu
      names(init) <- enum_dupl_string(rep_len("beta", p))
    } else {
      formula_parameters <- ~ts
      init <- rep_len(0, p * nsim)
      init[seq.int(from = 1L, by = nsim, length.out = p)] <- mu
      names(init) <- enum_dupl_string(rep_len("beta", p * nsim))
    }
  } else {
    formula_parameters <- ~(1 | ts)
    R <- chol(cov2cor(Sigma))
    iR <- upper.tri(R, diag = TRUE)
    R[iR] <- R[iR] * rep.int(1 / diag(R), seq_len(p))
    l <- list(
      beta = mu,
      theta = c(0.5 * log(diag(Sigma)), R[upper.tri(R, diag = FALSE)]),
      b = rep_len(0, nsim * p)
    )
    init <- unlist(l, FALSE, FALSE)
    names(init) <- enum_dupl_string(rep.int(names(l), lengths(l)))
  }
  environment(formula) <- environment(formula_windows) <-
    environment(formula_parameters) <- .GlobalEnv

  ## Create TMB object without optimizing
  zz <- egf(object,
    formula = formula,
    formula_windows = formula_windows,
    formula_parameters = formula_parameters,
    data = data,
    data_windows = data_windows,
    init = init,
    fit = FALSE
  )

  if (has_inflection && !is.null(Sigma)) {
    ## Second pass with 'data' of appropriate length,
    ## determined by simulated inflection times
    set_RNGstate()
    Y <- zz$tmb_out$simulate(init)$Y
    colnames(Y) <- names_top
    tmax <- ceiling(exp(Y[, "log(tinfl)"])) + 1
    time <- lapply(tmax, function(x) seq.int(0, x, by = 1))
    data <- data.frame(
      ts = rep.int(gl(nsim, 1L), lengths(time)),
      time = unlist(time, FALSE, FALSE),
      x = 0L # arbitrary non-negative integer to pass checks
    )
    data_windows <- data.frame(ts = gl(nsim, 1L), start = 0, end = tmax)
    zz <- update(zz, data = data, data_windows = data_windows)
  }

  ## Simulate
  set_RNGstate()
  sim <- zz$tmb_out$simulate(init)

  ## Replace dummy observations in 'data' with simulated ones
  data$x[] <- NA
  data$x[duplicated(data$ts)] <- sim$x

  ## Choose fitting window start times according to 'cstart' rule
  get_start <- function(d) {
    l <- c(0L, cumsum(d$x[-1L])) > cstart
    if (any(l)) d$time[which.max(l)] else NA_real_
  }
  data_windows$start <- c(by(data, data$ts, get_start))
  if (anyNA(data_windows$start)) {
    argna <- is.na(data_windows$start)
    data_windows$start[argna] <- 0
    warning(wrap(
      "Threshold 'cstart' not exceeded in these time series:\n\n",
      paste(sprintf("  %*d", as.integer(log10(nsim)) + 1L, which(argna)), collapse = "\n"), "\n\n",
      "Corresponding fitting windows contain all observations (for better or for worse)."
    ))
  }

  res <- list(
    model = object,
    mu = mu,
    Sigma = Sigma,
    Y = `colnames<-`(sim$Y, names_top),
    formula = formula,
    formula_windows = formula_windows,
    formula_parameters = formula_parameters,
    data = data,
    data_windows = data_windows,
    actual = init,
    nonrandom = grep("^(beta|theta)\\[", names(init)),
    call = match.call()
  )
  attr(res, "RNGstate") <- RNGstate
  class(res) <- c("egf_model_simulate", "list")
  res
}

#' Fit simulated incidence time series
#'
#' Estimates the generative model underlying a set of simulated
#' incidence time series.
#'
#' @param model
#'   An \code{"\link[=simulate.egf_model]{egf_model_simulate}"} object
#'   supplying simulated time series and details about the generative model
#'   to be estimated.
#' @param ...
#'   Optional arguments passed to \code{\link{egf.egf_model}},
#'   such as \code{priors} and \code{control}.
#'
#' @examples
#' example("simulate.egf_model", "epigrowthfit")
#' fit <- egf(sim)
#' pp <- data.frame(actual = sim$actual, fitted = fit$best)
#' pp[fit$nonrandom, , drop = FALSE]
#'
#' @export
egf.egf_model_simulate <- function(model, ...) {
  nel <- c("model", "formula", "formula_windows", "formula_parameters", "data", "data_windows")
  args <- model[nel]
  dots <- list(...)
  if (length(dots) > 0L && !is.null(nd <- names(dots))) {
    args <- c(args, dots[match(nd, nel, 0L) == 0L])
  }
  do.call(egf, args)
}
