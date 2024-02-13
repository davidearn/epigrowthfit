fitted.egf <-
function(object,
         top = egf_top(object),
         link = TRUE,
         se = FALSE,
         subset = NULL,
         select = NULL,
         ...) {
	stopifnot(is_true_or_false(link),
	          is_true_or_false(se))
	if (se && !link)
		stop("standard errors not available for inverse link-transformed fitted values")

	names_top <- egf_top(object)
	top <- unique(match.arg(top, names_top, several.ok = TRUE))

	frame_windows <- model.frame(object, "windows")
	frame_combined <- model.frame(object, "combined")
	subset <- egf_eval_subset(subset, frame_combined, parent.frame())
	select <- egf_eval_select(select, frame_combined, baseenv())

	if (se) {
		sdr <- egf_get_sdreport(object)
		ssdr <- summary(sdr, select = "report")
		index <- rownames(ssdr) == "Y"
		Y <- ssdr[index, "Estimate"]
		Y_se <- ssdr[index, "Std. Error"]
		dim(Y) <- dim(Y_se) <- object$tmb_out$env$ADreportDims$Y
	}
	else
		Y <- object$tmb_out$report(object$best)$Y

	## 'Y[i, j]' is the fitted value of top level nonlinear model parameter 'j'
	## (link scale) in fitting window 'i'
	colnames(Y) <- names_top
	Y <- Y[subset, top, drop = FALSE]

	res <- data.frame(top = rep(factor(top, levels = names_top),
	                            each = length(subset)),
	                  frame_windows[subset, c("ts", "window"), drop = FALSE],
	                  estimate = as.numeric(Y),
	                  row.names = NULL,
	                  check.names = FALSE)
	if (se) {
		colnames(Y_se) <- names_top
		Y_se <- Y_se[subset, top, drop = FALSE]
		res$se <- as.numeric(Y_se)
	}
	if (!link) {
		f <- lapply(egf_link_extract(levels(res$top)), egf_link_match,
		            inverse = TRUE)
		res$estimate <- in_place_ragged_apply(res$estimate, res$top, f = f)
		levels(res$top) <- egf_link_remove(levels(res$top))
	}
	res <- data.frame(res,
	                  frame_combined[subset, select, drop = FALSE],
	                  row.names = NULL,
	                  check.names = FALSE)
	attr(res, "se") <- se
	class(res) <- c("egf_fitted", oldClass(res))
	res
}

fitted.egf_no_fit <-
function(object,
         top = egf_top(object),
         link = TRUE,
         se = FALSE,
         subset = NULL,
         select = NULL,
         ...) {
	if (se)
		stop(gettextf("standard errors are not computed until model is estimated; retry with %s",
		              "object = update(object, se = TRUE, fit = TRUE, ...)"),
		     domain = NA)

	call <- match.call(expand.dots = FALSE)
	call[[1L]] <- quote(fitted.egf)
	object$best <- object$init
	eval(call)
}

confint.egf_fitted <-
function(object, parm, level = 0.95, link = TRUE, ...) {
	if (!attr(object, "se"))
		stop(gettextf("'%s' does not supply link scale fitted values and corresponding standard errors; retry with %s",
		              "object", "object = fitted(<egf>, link = TRUE, se = TRUE)"),
		     domain = NA)
	stopifnot(is_number_in_interval(level, 0, 1, "()"),
	          is_true_or_false(link))

	s <- c("top", "ts", "window", "estimate", "se")
	res <- data.frame(object[s[1:4]],
	                  wald(estimate = object$estimate, se = object$se,
	                       level = level),
	                  object[-match(s, names(object), 0L)],
	                  row.names = NULL,
	                  check.names = FALSE,
	                  stringsAsFactors = FALSE)
	if (!link) {
		elu <- c("estimate", "lower", "upper")
		f <- lapply(egf_link_extract(levels(res$top)), egf_link_match,
		            inverse = TRUE)
		res[elu] <- in_place_ragged_apply(res[elu], res$top, f = f)
		levels(res$top) <- egf_link_remove(levels(res$top))
	}
	attr(res, "level") <- level
	res
}
