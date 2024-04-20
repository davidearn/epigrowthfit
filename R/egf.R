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
	stopifnot(is.formula(formula_ts), is.formula(formula_windows))
	if (is.formula(formula_parameters))
		stopifnot(length(formula_parameters) == 2L)
	else
		stopifnot(is.list(formula_parameters),
		          all(vapply(formula_parameters, is.formula, FALSE)),
		          lengths(formula_parameters) == 3L)
	stopifnot(is.list(formula_priors),
	          all(vapply(formula_priors, is.formula, FALSE)),
	          lengths(formula_priors) == 3L)

	if (missing(data_ts))
		data_ts <- environment(formula_ts)
	else stopifnot(is.list(data_ts) || is.environment(data_ts))
	if (missing(data_windows))
		data_windows <- environment(formula_windows)
	else stopifnot(is.list(data_windows) || is.environment(data_windows))

	subset_ts <- substitute(subset_ts)
	subset_windows <- substitute(subset_windows)
	select_windows <- substitute(select_windows)

	na_action_ts <- match.arg(na_action_ts)
	na_action_windows <- match.arg(na_action_windows)

	stopifnot(inherits(control, "egf_control"),
	          isTrueFalse(fit), isTrueFalse(se),
	          is.list(init), is.list(map))

	names_parameters <- egf_top(model)

	formula_ts <- egf_sanitize_formula(formula_ts)
	formula_windows <- egf_sanitize_formula(formula_windows)
	formula_parameters <- egf_sanitize_formula_parameters(
		formula_parameters = formula_parameters,
		names_parameters = names_parameters,
		check_intercept = is.null(init))

	frame <- egf_make_frame(
		model = model,
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
	tmb_args <- egf_tmb_make_args(
		model = model,
		frame = frame,
		control = control,
		init = init,
		map = map,
		env = env)

	priors <- egf_make_priors(
		formula_priors = formula_priors,
		top   = list(names = names_parameters,
		             family = "norm"),
		beta  = list(length = env[["len"]][["beta"]],
		             family = "norm"),
		theta = list(length = env[["len"]][["theta"]],
		             family = "norm"),
		Sigma = list(length = length(tmb_args[["data"]][["block_rows"]]),
		             family = c("lkj", "wishart", "invwishart"),
		             rows = tmb_args[["data"]][["block_rows"]]))

	tmb_args[["data"]] <- egf_tmb_update_data(
		tmb_args[["data"]],
		priors = priors)

	tmb_out <- do.call(MakeADFun, tmb_args)
	tmb_out[["fn"]] <- egf_patch_fn(
		tmb_out[["fn"]],
		inner_optimizer = control[["inner_optimizer"]])
	tmb_out[["gr"]] <- egf_patch_gr(
		tmb_out[["gr"]],
		inner_optimizer = control[["inner_optimizer"]])

	ans <- list(model = model,
	            frame = frame,
	            priors = priors,
	            control = control,
	            tmb_out = tmb_out,
	            optimizer_out = NULL,
	            init = as.double(tmb_out[["env"]][["par"]]),
	            best = NULL,
	            random = tmb_out[["env"]][["lrandom"]]()
	            value = NULL,
	            gradient = NULL,
	            hessian = NULL,
	            effects = env[["effects"]],
	            contrasts = env[["contrasts"]],
	            call = match.call())

	if (!fit) {
		class(ans) <- "egf_no_fit"
		return(ans)
	}

	onomp <- openmp(n = NULL)
	if (on > 0L) {
		openmp(n = control[["omp_num_threads"]])
		on.exit(openmp(n = onomp))
	}
	optimizer <- control[["optimizer"]][["f"]]
	optimizer_args <- c(tmb_out[c("par", "fn", "gr")],
	                    control[["optimizer"]]["control"],
	                    control[["optimizer"]][["args"]])
	ans[["optimizer_out"]] <- do.call(optimizer, optimizer_args)

	ans[["best"]] <- as.double(tmb_out[["env"]][["last.par.best"]])
	ans[["value"]] <- as.double(tmb_out[["env"]][["value.best"]])
	if (se)
		rpt <- egf_adreport(ans, check = FALSE)
	if (!se || inherits(rpt, "error")) {
		ans[["gradient"]] <- tmb_out[["gr"]](ans[["best"]][!ans[["random"]]])
		ans[["hessian"]] <- NA
	}
	else {
		ans[["gradient"]] <- rpt[["gradient.fixed"]]
		ans[["hessian"]] <- rpt[["pdHess"]]
	}

	class(ans) <- "egf"
	ans
}
