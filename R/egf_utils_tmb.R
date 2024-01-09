##' Construct data objects for C++ template
##'
##' Gathers in a list data objects to be passed to
##' \code{\link[TMB]{MakeADFun}}.
##' These are retrieved inside of \pkg{epigrowthfit}'s C++ template
##' via \pkg{TMB}'s \code{DATA_*} macros.
##'
##' @param model
##'   An \code{"\link{egf_model}"} object.
##' @param frame
##'   A list returned by \code{egf_make_frame}.
##' @param control
##'   An \code{"\link{egf_control}"} object.
##' @param env
##'   An environment for storing intermediate objects that are not ultimately
##'   passed to \code{\link[TMB]{MakeADFun}}, but must nonetheless be preserved
##'   somewhere.
##'
##' @return
##' [Below,
##' \code{N = nlevels(frame$ts$window)}
##' is the number of fitting windows,
##' \code{n = N + sum(!is.na(frame$ts$window))}
##' is the total number of time points associated with a fitting window, and
##' \code{p = length(frame$parameters)}
##' is the number of top level nonlinear model parameters.]
##'
##' A list with elements:
##' \item{time}{
##'   A numeric vector of length \code{n} giving times since the left endpoint
##'   of the current fitting window.
##' }
##' \item{time_seg_len}{
##'   An integer vector of length \code{N} specifying the length of each
##'   fitting window as a number of time points.
##' }
##' \item{x}{
##'   A numeric vector of length \code{n-N} giving incidence in each
##'   fitting window. \code{x[i]} in window \code{k} is the number
##'   of cases observed from \code{time[k+i-1]} to \code{time[k+i]}.
##' }
##' \item{day1}{
##'   An integer vector of length \code{N}.
##'   If \code{model$day_of_week > 0}, then it indicates the first day
##'   of week in each fitting window, with value \code{i} in \code{0:6}
##'   mapping to the day of week \code{i} days after the reference day
##'   specified by \code{model$day_of_week}.
##'   Otherwise, it is filled with \code{-1}.
##' }
##' \item{flags}{
##'   A list with integer elements, used as flags to specify the model
##'   being estimated and to indicate what blocks of template code should
##'   be run.
##' }
##' \item{indices}{
##'   A list with integer elements and names of the form
##'   \code{"<link>_<parameter>"} (e.g., \code{"log_r"}),
##'   giving the column 0-index of top level nonlinear model parameters
##'   (e.g., \code{log(r)}) in the response matrix.
##'   Value \code{-1} is used for parameters not belonging to the model
##'   being estimated.
##' }
##' \item{Y}{
##'   The \link[=model.offset]{offset} component of the response matrix
##'   in dense format, with \code{N} rows and \code{p} columns.
##' }
##' \item{Xs, Xd}{
##'   The fixed effects design matrix in \link[Matrix:sparseMatrix]{sparse}
##'   or \link[=matrix]{dense} format, with \code{N} rows.
##'   If \code{control$sparse_X = TRUE}, then \code{Xs} is the sparse
##'   design matrix and \code{Xd} is an empty dense matrix.
##'   Otherwise, \code{Xs} is an empty sparse matrix and \code{Xd} is
##'   the dense design matrix.
##' }
##' \item{Z}{
##'   The random effects design matrix in \link[Matrix:sparseMatrix]{sparse}
##'   format, with \code{N} rows.
##'   If there are no random effects, then \code{Z} is an empty sparse matrix.
##' }
##' \item{beta_index, b_index}{
##'   Integer vectors of length \code{ncol(X)} and \code{ncol(Z)}, respectively,
##'   with values in \code{0:(p-1)}.
##'   These split the columns of \code{X} and \code{Z} by relation to a common
##'   top level nonlinear model parameter.
##' }
##' \item{beta_index_tab, b_index_tab}{
##'   Integer vectors of length \code{p} counting the columns of \code{X} and
##'   \code{Z}, respectively, that relate to a common top level nonlinear model
##'   parameter.
##' }
##' \item{block_rows, block_cols}{
##'   Integer vectors together giving the dimensions of each block of random
##'   effects coefficients.
##' }

egf_tmb_make_data <-
function(model, frame, control, env) {
	## Indices of time points associated with fitting windows
	first <- attr(frame$ts, "first")
	last <- attr(frame$ts, "last")
	index <- Map(seq.int, first, last)
	ulindex <- unlist1(index)

	## Fitting window lengths as numbers of time points
	time_seg_len <- lengths(index, use.names = FALSE)

	## Number of fitting windows
	N <- length(index)

	## Time since earliest time point
	time <- frame$ts$time[ulindex] - rep.int(frame$ts$time[first], time_seg_len)

	## Incidence
	x <- frame$ts$x[!is.na(frame$ts$window)]

	if (model$day_of_week > 0L) {
		## Date during 24 hours starting at earliest time point
		Date1 <- .Date(frame$ts$time[first])

		## Day of week on that Date coded as an integer 'i' in '0:6'.
		## Integer 'i' maps to the day of week 'i' days after a reference day,
		## which is the day 'model$day_of_week' days after Saturday.
		## Less verbosely: i -> weekdays(.Date(2) + model$day_of_week + i),
		## noting that weekdays(.Date(2)) == "Saturday" ...
		origin <- .Date(2L + model$day_of_week)
		day1 <- as.integer(julian(Date1, origin = origin) %% 7)
	} else {
		day1 <- rep.int(-1L, N)
	}

	## Response matrix, offset component only
	offsets <- lapply(frame$parameters, model.offset)
	offsets[vapply(offsets, is.null, FALSE)] <- list(double(N))
	Y <- do.call(cbind, offsets)

	## Names of top level nonlinear model parameters
	names_top <- names(frame$parameters)

	## List of fixed effects formulae and list of lists of random effects terms
	l <- lapply(frame$parameters, function(x) split_effects(formula(terms(x))))
	fixed <- lapply(l, `[[`, "fixed")
	random <- lapply(l, `[[`, "random")

	## Fixed effects infrastructure
	X <- Map(egf_make_X, fixed = fixed,
	         data = frame$parameters, sparse = control$sparse_X)
	Xc <- egf_combine_X(fixed = fixed, X = X)
	beta_index <- as.integer(Xc$effects$top) - 1L
	beta_index_tab <- c(table(Xc$effects$top))

	## Random effects infrastructure
	random1 <- unlist1(random)
	names(random1) <- rep.int(names(random), lengths(random))
	Z <- Map(egf_make_Z, random = random1,
	         data = rep.int(frame$parameters, lengths(random)))
	Zc <- egf_combine_Z(random = random1, Z = Z)
	if (is.null(Zc)) {
		b_index <- integer(0L)
		b_index_tab <- integer(length(names_top))
		block_rows <- integer(0L)
		block_cols <- integer(0L)
	} else {
		Zc$effects$top <- factor(Zc$effects$top, levels = names_top)
		b_index <- as.integer(Zc$effects$top) - 1L
		b_index_tab <- c(table(Zc$effects$top))
		block_rows <-
			as.integer(colSums(table(Zc$effects$top, Zc$effects$cov) > 0L))
		block_cols <-
			as.integer(colSums(table(Zc$effects$vec, Zc$effects$cov) > 0L))
	}

	## Flags
	flags <- list(curve = egf_get_flag(model$curve, "curve"),
	              excess = as.integer(model$excess),
	              family = egf_get_flag(model$family, "family"),
	              day_of_week = as.integer(model$day_of_week > 0L),
	              trace = control$trace,
	              sparse_X = as.integer(control$sparse_X),
	              predict = 0L)

	## Column indices of top level nonlinear model parameters
	## in response matrix
	names_top_all <- egf_get_names_top(NULL, link = TRUE)
	indices <- as.list(match(names_top_all, names_top, 0L) - 1L)
	names(indices) <- sub("^(log|logit)\\((.*)\\)$", "\\1_\\2", names_top_all)

	## Stuff to preserve but not return
	env$effects <- list(beta = Xc$effects, b = Zc$effects)
	env$contrasts <- list(X = Xc$contrasts, Z = Zc$contrasts)
	env$len <- c(beta = ncol(Xc$X),
	             theta = sum(as.integer(choose(block_rows + 1L, 2L))),
	             b = sum(block_rows * block_cols))

	list(time = time,
	     time_seg_len = time_seg_len,
	     x = x,
	     day1 = day1,
	     flags = flags,
	     indices = indices,
	     Y = Y,
	     Xd = if (control$sparse_X) matrix(double(0L), N, 0L)
	          else Xc$X,
	     Xs = if (control$sparse_X) Xc$X
	          else Matrix(double(0L), N, 0L, sparse = TRUE),
	     Z = if (is.null(Zc)) Matrix(double(0L), N, 0L, sparse = TRUE)
	         else Zc$Z,
	     beta_index = beta_index,
	     b_index = b_index,
	     beta_index_tab = beta_index_tab,
	     b_index_tab = b_index_tab,
	     block_rows = block_rows,
	     block_cols = block_cols)
}

##' Construct parameter objects for C++ template
##'
##' Gathers in a list parameter objects to be passed to
##' \code{\link[TMB]{MakeADFun}}
##' and used during the first likelihood evaluation.
##' These are retrieved inside of \pkg{epigrowthfit}'s C++ template
##' via \pkg{TMB}'s \code{PARAMETER_*} macros.
##'
##' @param model
##'   An \code{"\link{egf_model}"} object.
##' @param frame
##'   A list returned by \code{egf_make_frame}.
##' @param env
##'   An environment for storing intermediate objects that are not ultimately
##'   passed to \code{\link[TMB]{MakeADFun}}, but must nonetheless be preserved
##'   somewhere.
##'
##' @details
##' Naive estimates of top level nonlinear model parameters are obtained
##' for each fitting window as follows:
##' \describe{
##' \item{\code{r}}{
##'   The slope of a linear model fit to \code{log1p(cumsum(x)))}.
##' }
##' \item{\code{c0}}{
##'   \code{exp(log_c0)}, where \code{log_c0} is the intercept of
##'   a linear model fit to \code{log1p(cumsum(x))}.
##' }
##' \item{\code{tinfl}}{
##'   \code{max(time)}. This assumes that the fitting window ends
##'   near the time of a peak in incidence.
##' }
##' \item{\code{K}}{
##'   \code{2*sum(x)}. This assumes that the fitting window ends
##'   near the time of a peak in incidence _and_ that incidence
##'   is roughly symmetric about the peak.
##' }
##' \item{\code{p}}{0.95}
##' \item{\code{a, b, disp, w[123456]}}{1}
##' \item{\code{alpha}}{
##'   \code{r*c0^(1-p)} if \code{curve = "subexponential"},
##'   \code{r/log(K/c0)} if \code{curve = "gompertz"}.
##'   These are the values obtained by setting the per capita growth
##'   rate at time 0 in the subexponential and Gompertz models equal
##'   to \code{r}, substituting the naive estimates of \code{r},
##'   \code{c0}, \code{K}, and \code{p}, and solving for \code{alpha}.
##' }
##' }
##' The naive estimates are link-transformed, and the means of the
##' link scale estimates across fitting windows are used as initial
##' values for corresponding \code{"(Intercept)"} coefficients in
##' \code{beta}.
##'
##' @return
##' A list with elements \code{beta}, \code{theta}, and \code{b},
##' each numeric vectors. \code{theta} and \code{b} are zero vectors,
##' while \code{beta} is a zero vector except for \code{"(Intercept)"}
##' coefficients; see Details.

egf_tmb_make_parameters <-
function(model, frame, env) {
	## Initialize each parameter object to a zero vector
	res <- lapply(env$len, double)

	## Identify top level nonlinear model parameters
	## whose mixed effects formula has an intercept
	f <- function(x) attr(terms(x), "intercept") == 1L
	has1 <- vapply(frame$parameters, f, FALSE)
	if (!any(has1)) {
		return(res)
	}

	## Record naive estimates of top level nonlinear model parameters
	## (one per time series segment)
	tx <- split(frame$ts[c("time", "x")], frame$ts$window)
	compute_naive <-
    function(d) {
		n <- max(2, trunc(nrow(d) / 2))
		ab <- try(coef(lm(log1p(cumsum(x)) ~ time,
		                  data = d,
		                  subset = seq_len(n),
		                  na.action = na.omit)),
		          silent = TRUE)
		if (inherits(ab, "try-error") || !all(is.finite(ab))) {
			r <- 0.04
			c0 <- 1
		} else {
			r <- ab[[2L]]
			c0 <- exp(ab[[1L]])
		}
		tinfl <- max(d$time)
		K <- 2 * sum(d$x, na.rm = TRUE)
		c(r = r, c0 = c0, tinfl = tinfl, K = K)
	}
	naive <- as.data.frame(t(vapply(tx, compute_naive, double(4L))))
	naive["p"] <- list(0.95)
	naive[c("a", "b", "disp", paste0("w", 1:6))] <- list(1)
	naive$alpha <- switch(model$curve,
	                      subexponential = naive$r * naive$c0^(1 - naive$p),
	                      gompertz = naive$r / (log(naive$K) - log(naive$c0)))

	## Link transform
	link <- lapply(egf_link_get(names(naive)), egf_link_match)
	naive[] <- Map(function(f, x) f(x), link, naive)
	names(naive) <- egf_link_add(names(naive))

	## Identify elements of 'beta' corresponding to "(Intercept)"
	## and assign means of naive estimates over time series segments
	names_top <- names(frame$parameters)
	m <- match(names_top[has1], env$effects$beta$top, 0L)
	res$beta[m] <- colMeans(naive[names_top[has1]])
	res
}

egf_tmb_make_args <-
function(model, frame, control, init, map, env) {
	data <- egf_tmb_make_data(model = model, frame = frame, control = control,
	                          env = env)
	parameters <- egf_tmb_make_parameters(model = model, frame = frame,
	                                      env = env)

	nms <- names(parameters)
	init <- init[match(nms, names(init), 0L)]
	map <- map[match(nms, names(map), 0L)]
	for (s in nms) {
		n <- length(parameters[[s]])
		if (s %in% names(init)) {
			eval(bquote(stopifnot(is.numeric(init[[.(s)]]),
			                      length(init[[.(s)]]) == .(n),
			                      !is.infinite(init[[.(s)]]))))
			## Replace
			index <- !is.na(init[[s]])
			parameters[[s]][index] <- init[[s]][index]
		}
		if (s %in% names(map)) {
			if (is.factor(map[[s]])) {
				eval(bquote(stopifnot(length(map[[.(s)]]) == .(n))))
				map[[s]] <- factor(map[[s]])
			} else {
				eval(bquote(stopifnot(!anyNA(index <- seq_len(.(n))[map[[.(s)]]]))))
				f <- rep.int(NA_integer_, n)
				f[-index] <- seq_len(n - length(index))
				levels(f) <- as.character(f[-index])
				class(f) <- "factor"
				map[[s]] <- f
			}
		}
	}

	if (length(parameters$b) > 0L) {
		## Declare that 'b' contains random effects
		random <- "b"
	} else {
		## Declare that there are no random effects
		random <- NULL
		## Fix 'theta' and 'b' to NA_real_ since only 'beta' is used
		## and zero-length vectors are disallowed
		parameters$theta <- parameters$b <- NA_real_
		map$theta <- map$b <- factor(NA)
	}

	list(data = data,
	     parameters = parameters,
	     map = map,
	     random = random,
	     profile = if (control$profile) "beta" else NULL,
	     DLL = "epigrowthfit",
	     silent = (control$trace == 0L))
}

egf_tmb_update_data <-
function(data, priors) {
	priors$bottom <- unlist1(priors$bottom)
	for (s in c("top", "bottom")) {
		l <- priors[[s]]
		n <- length(l)
		argnull <- vapply(l, is.null, FALSE)

		flag <- egf_get_flag(vapply(l[!argnull], `[[`, "", "family"), "prior")
		flag <- replace(rep.int(-1L, n), !argnull, flag)
		data$flags[[paste0("regularize_", s)]] <- flag

		par <- lapply(l[!argnull], function(x) unlist1(x$parameters))
		par <- replace(rep.int(list(double(0L)), n), !argnull, par)
		data[[paste0("hyperparameters_", s)]] <- par
	}
	data
}

egf_tmb_remake_args <-
function(obj, par) {
	nms <- c("data", "parameters", "map", "random", "profile", "DLL", "silent")
	res <- mget(nms, envir = obj$env)
	res$parameters <- obj$env$parList(par[obj$env$lfixed()], par)
	attr(res$data, "check.passed") <- NULL
	attr(res$parameters, "check.passed") <- NULL
	if (ncol(res$data$Z) > 0L) {
		res$random <- "b"
	}
	if (!is.null(res$profile)) {
		res$profile <- "beta"
	}
	res
}
