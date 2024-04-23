## MJ: priority for refactoring

## NOTE: This evalates cumulative incidence as a function of time
##       conditional on an estimated nonlinear mixed effects model.
##       Originally, it was registered as a method for generic function
##       'predict', but that was unconventional and now it is an
##       internal function used only by the 'plot' method.  Hence:
##
## TODO: Implement a _proper_ method for 'predict' and figure out
##       how best to expose whatever _this_ is ...

.__predict.egf <-
function(object,
         what = c("interval", "cumulative", "rt"),
         time,
         window,
         log = TRUE,
         se = FALSE,
         ...) {
	what <- unique(match.arg(what, several.ok = TRUE))
	stopifnot(isTrueFalse(log), isTrueFalse(se))
	if (se && !log)
		stop("standard errors not available for inverse log-transformed predicted values")

	frame_windows <- model.frame(object, "windows")
	start <- frame_windows$start
	end <- frame_windows$end
	day1 <- object$tmb_out$env$data$day1
	do_day_of_week <- object$model$day_of_week > 0L

	if (missing_time <- missing(time)) {
		len  <- object$tmb_out$env$data$time_seg_len
		subset <- seq_along(len)
		time <- object$tmb_out$env$data$time + rep.int(start, len)
		time_split <- split(time, rep.int(subset, len))
	}
	else {
		if (inherits(time, c("Date", "POSIXct", "POSIXlt")))
			time <- julian(time)
		stopifnot(is.numeric(time), length(time) > 0L, is.finite(time))
		if (do_day_of_week) {
			stopifnot(all.equal(time, z <- round(time)))
			time <- z
		}
		stopifnot(is.factor(window), length(window) == length(time))

		subset <- which(levels(frame_windows$window)%in%levels(factor(window)))
		window <- factor(window, levels = levels(frame_windows$window)[subset])

		if (nlevels(window) == 0L)
			stop(gettextf("'%s' must have at least one valid level", "window"),
			     domain = NA)
		len <- c(table(window))
		min_len <- 1L + as.integer("interval" %in% what)
		if (any(len < min_len))
			stop(gettextf("'%s' does not have minimum length %d in each level of '%s'",
			              "time", min_len, "window"),
			     domain = NA)

		time_split <- split(time, window)
		t0 <- vapply(time_split, min, 0)
		t1 <- vapply(time_split, max, 0)

		start <- start[subset]
		end <- end[subset]
		day1 <- day1[subset]

		if (any(t0 < start | t1 > end))
			stop(gettextf("%s[i] occurs before the start or after the end of %s[i]",
			              "time", "window"),
			     domain = NA)
		check_ok_diff_time <-
			if (do_day_of_week)
				function(x) all(diff(x) == 1)
			else
				function(x) all(diff(x) > 0)
		if (!all(vapply(time_split, check_ok_diff_time, FALSE)))
			stop(switch(1L + do_day_of_week,
			            gettextf("'%s' must be increasing in each level of '%'",
			                     "time", "window"),
			            gettextf("'%s' must be increasing with unit spacing in each level of '%'",
			                     "time", "window")),
			     domain = NA)

		time <- unlist1(time_split)
		if (do_day_of_week)
			day1 <- as.integer((day1 + (t0 - start)) %% 7)
	}

	tmb_args <- egf_tmb_remake_args(object$tmb_out, par = object$best)
	tmb_args$data$flags$predict <- 1L
	tmb_args$data$what <- as.integer(eval(formals(sys.function())$what)%in%what)
	tmb_args$data$subset <- subset - 1L
	tmb_args$data$new_time <- time - rep.int(start, len)
	tmb_args$data$new_time_seg_len <- len
	tmb_args$data$new_day1 <- day1
	tmb_out_retape <- do.call(MakeADFun, tmb_args)

	if (se) {
		sdr <- sdreport(tmb_out_retape,
		                par.fixed = object$best[!object$random],
		                getReportCovariance = FALSE)
		ssdr <- summary(sdr, select = "report")
		index <- factor(rownames(ssdr),
		                levels = sprintf("log_%s", what),
		                labels = sprintf("log(%s)", what))
		report <- split(unname(ssdr[, "Estimate"]), index)
		report_se <- split(unname(ssdr[, "Std. Error"]), index)
	}
	else {
		report <- tmb_out_retape$report(object$best)[sprintf("log_%s", what)]
		names(report) <- sprintf("log(%s)", what)
	}

	last <- cumsum(len)
	first <- c(0L, last[-length(last)]) + 1L
	x <- rep.int(NA_real_, length(time))
	ix <- list(interval = -first,
	           cumulative = seq_along(time),
	           rt = seq_along(time))

	res <- data.frame(var = gl(length(report), length(time),
	                           labels = names(report)),
	                  ts = rep.int(frame_windows$ts[subset], len),
	                  window = rep.int(frame_windows$window[subset], len),
	                  time = time,
	                  value =
	                      unlist1(Map(replace, list(x), ix[what], report)))
	if (se)
		res$se <- unlist1(Map(replace, list(x), ix[what], report_se))
	if (!log) {
		res$value <- exp(res$value)
		levels(res$var) <- what
	}
	attr(res, "se") <- se
	class(res) <- c("predict.egf", oldClass(res))
	res
}

.__confint.predict.egf <-
function(object, parm, level = 0.95, log = TRUE, ...) {
	stopifnot(isNumber(level), level > 0, level < 1,
	          isTrueFalse(log), attr(object, "se"))

	res <- data.frame(object[-match("se", names(object), 0L)],
	                  wald(object$value, object$se, level = level))
	attr(res, "level") <- level
	if (!log) {
		elu <- c("value", "lower", "upper")
		res[elu] <- exp(res[elu])
		levels(res$var) <- sub("^log\\((.+)\\)$", "\\1", levels(res$var))
	}
	res
}
