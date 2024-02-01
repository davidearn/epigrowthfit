egf <-
function(model, ...)
	UseMethod("egf", model)

egf.egf_model <-
function(model,
         formula_ts,
         formula_windows,
         formula_parameters = list(),
         formula_priors = list(),
         data_ts,
         data_windows,
         subset_ts = NULL,
         subset_windows = NULL,
         select_windows = NULL,
         na_action_ts = c("fail", "pass"),
         na_action_windows = c("fail", "omit"),
         control = egf_control(),
         fit = TRUE,
         se = FALSE,
         init = list(),
         map = list(),
         ...) {
	stopifnot(inherits(formula_ts, "formula"),
	          inherits(formula_windows, "formula"))
	if (inherits(formula_parameters, "formula")) {
		stopifnot(length(formula_parameters) == 2L)
	} else {
		stopifnot(is.list(formula_parameters),
		          vapply(formula_parameters, inherits, FALSE, "formula"),
		          lengths(formula_parameters) == 3L)
	}
	stopifnot(is.list(formula_priors),
	          vapply(formula_priors, inherits, FALSE, "formula"),
	          lengths(formula_priors) == 3L)
	if (missing(data_ts)) {
		data_ts <- environment(formula_ts)
	} else {
		stopifnot(is.list(data_ts) || is.environment(data_ts))
	}
	if (missing(data_windows)) {
		data_windows <- environment(formula_windows)
	} else {
		stopifnot(is.list(data_windows) || is.environment(data_windows))
	}
	subset_ts <- substitute(subset_ts)
	subset_windows <- substitute(subset_windows)
	select_windows <- substitute(select_windows)
	na_action_ts <- match.arg(na_action_ts)
	na_action_windows <- match.arg(na_action_windows)
	stopifnot(inherits(control, "egf_control"),
	          is_true_or_false(fit),
	          is_true_or_false(se),
	          is.list(init),
	          is.list(map))

	names_parameters <- egf_get_names_top(model, link = TRUE)

	formula_ts <- egf_sanitize_formula(formula_ts)
	formula_windows <- egf_sanitize_formula(formula_windows)
	formula_parameters <-
		egf_sanitize_formula_parameters(formula_parameters = formula_parameters,
		                                names_parameters = names_parameters,
		                                check_intercept = is.null(init))
	frame <- egf_make_frame(model = model,
	                        formula_ts = formula_ts,
	                        formula_windows = formula_windows,
	                        formula_parameters = formula_parameters,
	                        data_ts = data_ts,
	                        data_windows = data_windows,
	                        subset_ts = subset_ts,
	                        subset_windows = subset_windows,
	                        select_windows = select_windows,
	                        na_action_ts = na_action_ts,
	                        na_action_windows = na_action_windows)

	env <- new.env(parent = emptyenv())
	tmb_args <- egf_tmb_make_args(model = model,
	                              frame = frame,
	                              control = control,
	                              init = init,
	                              map = map,
	                              env = env)

	priors <-
		egf_make_priors(formula_priors = formula_priors,
		                top = list(names = names_parameters,
		                           family = "norm"),
		                beta = list(length = env$len[["beta"]],
		                            family = "norm"),
		                theta = list(length = env$len[["theta"]],
		                             family = "norm"),
		                Sigma = list(length = length(tmb_args$data$block_rows),
		                             family = c("lkj", "wishart", "invwishart"),
		                             rows = tmb_args$data$block_rows))
	tmb_args$data <- egf_tmb_update_data(tmb_args$data, priors = priors)

	tmb_out <- do.call(MakeADFun, tmb_args)
	tmb_out$fn <- egf_patch_fn(tmb_out$fn,
	                           inner_optimizer = control$inner_optimizer)
	tmb_out$gr <- egf_patch_gr(tmb_out$gr,
	                           inner_optimizer = control$inner_optimizer)

	res <- list(model = model,
	            frame = frame,
	            priors = priors,
	            control = control,
	            tmb_out = tmb_out,
	            optimizer_out = NULL,
	            init = tmb_out$env$par,
	            best = NULL,
	            random = tmb_out$env$lrandom(),
	            value = NULL,
	            gradient = NULL,
	            hessian = NULL,
	            sdreport = NULL,
	            effects = env$effects,
	            contrasts = env$contrasts,
	            call = match.call())

	if (!fit) {
		class(res) <- "egf_no_fit"
		return(res)
	}

	on <- openmp(n = NULL)
	if (on > 0L) {
		openmp(n = control$omp_num_threads)
		on.exit(openmp(n = on))
	}
	optimizer <- control$optimizer$f
	optimizer_args <- c(tmb_out[c("par", "fn", "gr")],
	                    control$optimizer["control"],
	                    control$optimizer[["args"]])
	res$optimizer_out <- do.call(optimizer, optimizer_args)

	res$best <- tmb_out$env$last.par.best
	res$value <- as.double(tmb_out$env$value.best)
	if (se) {
		res$sdreport <- try(sdreport(tmb_out,
		                             par.fixed = res$best[!res$random],
		                             getReportCovariance = FALSE))
	}
	if (inherits(res$sdreport, "sdreport")) {
		res$gradient <- res$sdreport$gradient.fixed
		res$hessian <- res$sdreport$pdHess
	} else {
		res$gradient <- tmb_out$gr(res$best[!res$random])
		res$hessian <- NA
	}

	class(res) <- "egf"
	res
}
