## MJ: priority for refactoring

plot.egf <-
function(x,
         type = c("interval", "cumulative", "rt"),
         time_as = c("Date", "numeric"),
         delta = 1,
         log = TRUE,
         zero = NA,
         show_predict = TRUE,
         show_doubling = FALSE,
         level = 0.95,
         control = egf_control_plot(),
         cache = NULL,
         plot = TRUE,
         subset = NULL,
         order = NULL,
         xlim = NULL,
         ylim = NULL,
         main = NULL,
         sub = NULL,
         xlab = NULL,
         ylab = NULL,
         ...) {
	## Validate arguments ------------------------------------------------------

	type <- match.arg(type)
	if (x$model$day_of_week > 0L)
		delta <- 1
	else stopifnot(isNumber(delta), delta > 0)
	stopifnot(isInteger(show_predict))
	show_predict <- min(2L, max(0L, as.integer(show_predict)))
	if (x$model$curve %in% c("exponential", "logistic", "richards")) {
		stopifnot(isInteger(show_doubling))
		show_doubling <- min(2L, max(0L, as.integer(show_doubling)))
		show_asymptote <-
			as.integer(type == "rt" && x$model$curve != "exponential")
	}
	else {
		show_doubling <- 0L
		show_asymptote <- 0L
	}
	stopifnot(isTrueFalse(plot))
	if (plot) {
		time_as <- match.arg(time_as)
		stopifnot(isTrueFalse(log))
		if (!is.null(zero)) {
			if (type == "rt")
				zero <- NULL
			else if (is.na(zero))
				zero <- as.double(zero)
			else stopifnot(isNumber(zero), zero > 0)
		}
		if (any(c(show_predict, show_doubling) == 2L))
			stopifnot(isNumber(level), level > 0, level < 1)
		stopifnot(inherits(control, "egf_control_plot"))
		if (!is.null(cache))
			stopifnot(inherits(cache, "plot.egf"))
		if (!is.null(xlim)) {
			if (!is.numeric(xlim))
				xlim <- try(julian(as.Date(xlim)))
			stopifnot(is.numeric(xlim),
			          length(xlim) == 2L,
			          is.finite(xlim),
			          xlim[1L] < xlim[2L])
		}
		if (!is.null(ylim))
			stopifnot(is.numeric(ylim),
			          length(ylim) == 2L,
			          is.finite(ylim),
			          ylim[1L] < ylim[2L])
	}


	## Subset and order time series, fitting windows ---------------------------

	frame_ts <- model.frame(x, full = TRUE)
	frame_windows <- frame_windows_bak <- model.frame(x, "windows")
	frame_combined <- model.frame(x, "combined")

	subset <- egf_eval_subset(subset, frame_combined, parent.frame())
	if (length(subset) == 0L)
		stop(gettextf("'%s' is empty; nothing to plot", "subset"),
		     domain = NA)
	order <- egf_eval_order(order, frame_combined, parent.frame())
	subset <- order[order %in% subset]
	frame_windows <- frame_windows[subset, , drop = FALSE]

	lts <- as.character(unique(frame_windows$ts))
	frame_windows$ts <- factor(frame_windows$ts, levels = lts)

	reorder <- do.call(base::order, unname(frame_windows[c("ts", "start")]))
	subset <- subset[reorder]
	frame_windows <- frame_windows[reorder, , drop = FALSE]

	lw <- as.character(frame_windows$window)
	frame_windows$window <- factor(frame_windows$window, levels = lw)

	frame_ts$ts <- factor(frame_ts$ts, levels = lts)
	frame_ts$window <- factor(frame_ts$window, levels = lw)
	frame_ts <- frame_ts[!is.na(frame_ts$ts), , drop = FALSE]


	## Finalize annotation -----------------------------------------------------

	if (plot) {
		n <- nlevels(frame_ts$ts)
		subset1 <- subset[match(lts, frame_windows$ts, 0L)]
		nplot <- n
		main <- egf_eval_label(main, frame_combined, parent.frame())[subset1]
		if (is.null(main))
			main <- levels(frame_ts$ts)
		sub <- egf_eval_label(sub, frame_combined, parent.frame())[subset1]
		if (is.null(sub))
			sub <- ""
		main <- rep_len(main, nplot)
		sub <- rep_len(sub, nplot)

		if (is.null(xlab))
			xlab <- switch(time_as, Date = "", numeric = "time")
		if (type == "rt") {
			if (is.null(ylab))
				ylab <- c("per capita growth rate", "doubling time")
			else if ((is.character(ylab) || is.expression(ylab)) &&
			         length(ylab) == 1L)
				ylab <- c(ylab, "doubling time")
		}
		else {
			if (is.null(ylab))
				ylab <- paste(type, "incidence")
		}
	}


	## Augment 'cache' ---------------------------------------------------------

	## If necessary, initialize
	if (is.null(cache))
		cache <- data.frame(var = factor(),
		                    ts = factor(),
		                    window = factor(),
		                    time = double(0L),
		                    value = double(0L),
		                    se = double(0L))
	nc <- names(cache)

	## If necessary, augment 'cache' with predicted values
	## of whatever is being plotted and standard errors
	if (show_predict > 0L) {
		ok <-
			cache$var == sprintf("log(%s)", type) &
			!(show_predict == 2L & is.na(cache$se))
		required <- setdiff(lw, cache$window[ok])
		if (length(required) > 0L) {
			m <- match(required, frame_windows_bak$window, 0L)
			time_split <- Map(seq.int,
			                  from = frame_windows_bak$start[m],
			                  to = frame_windows_bak$end[m],
			                  by = delta)
			px <- .__predict.egf(x,
			              what = type,
			              time = unlist1(time_split),
			              window = rep.int(frame_windows_bak$window[m],
			                               lengths(time_split)),
			              log = TRUE,
			              se = (show_predict == 2L))
			if (show_predict == 1L)
				px$se <- NA_real_
			cache <- rbind(cache, px[nc])
		}
	}

	## If necessary, augment 'cache' with fitted values of 'log(r)'
	## and standard errors
	if (show_doubling > 0L || show_asymptote > 0L) {
		ok <- cache$var == "log(r)" & !(show_doubling == 2L & is.na(cache$se))
		required <- setdiff(lw, cache$window[ok])
		if (length(required) > 0L) {
			fx <- fitted(x,
			             top = "log(r)",
			             class = TRUE,
			             se = (show_doubling == 2L),
			             subset = (frame_windows_bak$window %in% required))
			if (show_doubling != 2L)
				fx$se <- NA_real_
			fx$time <- NA_real_
			names(fx)[match("top", names(fx), 0L)] <- "var"
			cache <- rbind(cache, fx[nc])
		}
	}

	## Clean up
	o <- do.call(base::order, unname(cache[-5L]))
	cache <- cache[o, , drop = FALSE]
	i <- !duplicated(cache[1:4])
	cache <- cache[i, , drop = FALSE]
	row.names(cache) <- NULL
	class(cache) <- c("plot.egf", oldClass(cache))

	## If not plotting, then return
	if (!plot)
		return(invisible(cache))

	## If plotting, then create an instruction to return
	## 'cache' if the low level plot function throws an error
	cache_bak <- cache
	on.exit({
		message(gettextf("augmented '%s' returned despite error", "cache"),
		        domain = NA)
		return(invisible(cache_bak))
	})

	## Extract only those rows of 'cache' needed by the low level plot function
	lvar <- c(sprintf("log(%s)", type), if (show_doubling > 0L) "log(r)")
	cache[1:3] <- Map(factor, cache[1:3], levels = list(lvar, lts, lw))
	i <- complete.cases(cache[1:3])
	cache <- cache[i, , drop = FALSE]

	## Augment with confidence intervals
	cache[c("lower", "upper")] <- list(NA_real_)
	i <-
		(show_predict == 2L & cache$var == sprintf("log(%s)", type)) |
		(show_doubling == 2L & cache$var == "log(r)")
	if (any(i))
		cache[i, c("lower", "upper")] <-
			do.call(wald, c(cache[i, c("value", "se"), drop = FALSE],
			                list(level = level)))


	## Augment 'control' -------------------------------------------------------

	## Gather parameter values necessary to determine text dimensions ...
	## those not specified in 'control' are obtained from 'gp'
	par(cex = 1)
	gp <- par()
	if (is.list(control$axis$y)) {
		nms <- setdiff(c("mgp", "cex.axis", "font.axis", "family"),
		               names(control$axis$y))
		control$axis$y[nms] <- gp[nms]
	}
	if (is.list(control$title$sub)) {
		nms <- setdiff(c("mgp", "cex.sub", "font.sub", "family"),
		               names(control$title$sub))
		control$title$sub[nms] <- gp[nms]
	}
	if (is.list(control$title$ylab)) {
		nms <- setdiff(c("mgp", "cex.lab", "font.lab", "family"),
		               names(control$title$ylab))
		control$title$ylab[nms] <- gp[nms]
	}
	for (s in c("ci", "value", "legend"))
		if (is.list(control$doubling[[s]])) {
			nms <- setdiff(c("adj", "cex", "font", "family"),
			               names(control$doubling[[s]]))
			control$doubling[[s]][nms] <- gp[nms]
		}

	## Distinguish minor and major axes ...
	if (time_as == "Date" && is.list(args <- control$axis$x)) {
		nms <- setdiff(c("mgp", "cex.axis"), names(args))
		args[nms] <- gp[nms]
		args <- rep.int(list(args), 2L)
		names(args) <- c("minor", "major")
		args$major$mgp[1:2] <- args$minor$mgp[1:2] + 1.5
		args$major$cex.axis <- args$minor$cex.axis * 1.15
		args$major$tick <- FALSE
		control$axis$x <- args
	}


	## Misc --------------------------------------------------------------------

	## Approximation of per capita growth rate from cumulative counts
	## is much less noisy but sensitive to the (unknown) initial value,
	## so retrieve the predicted initial value and hope that it's reasonable
	if (type == "rt") {
		px <- .__predict.egf(x,
		              what = "cumulative",
		              time = frame_windows$start,
		              window = frame_windows$window,
		              log = FALSE,
		              se = FALSE)
		frame_windows$c0 <-
			px$value[match(frame_windows$window, px$window, 0L)]
	}


	## Plot and return ---------------------------------------------------------

	plot.egf.curve(frame_ts = frame_ts,
	               frame_windows = frame_windows,
	               cache = cache,
	               type = type,
	               time_as = time_as,
	               delta = delta,
	               log = log,
	               zero = zero,
	               show_predict = show_predict,
	               show_doubling = show_doubling,
	               show_asymptote = show_asymptote,
	               level = level,
	               control = control,
	               xlim = xlim,
	               ylim = ylim,
	               main = main,
	               sub = sub,
	               xlab = xlab,
	               ylab = ylab)

	## Discard exit instructions if low level plot function runs without error
	on.exit()
	invisible(cache)
}

plot.egf.curve <-
function(frame_ts, frame_windows, cache,
         type, time_as, delta, log, zero,
         show_predict, show_doubling, show_asymptote,
         level, control, xlim, ylim, main, sub, xlab, ylab) {
	n <- nlevels(frame_ts$ts)
	formula <- as.formula(call("~", as.name(type), quote(time)))
	xlim_bak <- xlim
	ylim_bak <- ylim
	elu <- c("value", "lower", "upper")

	## Graphical parameters
	gp <- par()

	for (i in seq_len(n)) { # loop over plots
		data <- frame_ts[unclass(frame_ts$ts) == i, , drop = FALSE]
		data_windows <- frame_windows[unclass(frame_windows$ts) == i, , drop = FALSE]
		if (show_doubling > 0L) {
			cache_r <- cache[unclass(cache$ts) == i & cache$var == "log(r)", , drop = FALSE]
			cache_r[elu] <- exp(cache_r[elu])
		}
		if (show_predict > 0L) {
			cache_predict <- cache[unclass(cache$ts) == i & cache$var == sprintf("log(%s)", type), , drop = FALSE]
			cache_predict[elu] <- exp(cache_predict[elu])
		}

		data$delta <- c(NA, diff(data$time))
		data[[type]] <-
			switch(type,
			       interval = c(NA, data$x[-1L]) * delta / data$delta,
			       cumulative = c(0, cumsum(data$x[-1L])),
			       rt = local({
			       	f <-
			       	function(x, x0)
			       		c(0.5 * diff(log(cumsum(c(x0, x))), lag = 2), NA)
			       	y <- unsplit(Map(f,
			       	                 x = split(data$x,
			       	                           data$window,
			       	                           drop = TRUE),
			       	                 x0 = data_windows$c0),
			       	             data$window,
			       	             drop = TRUE)
			       	y[!is.finite(y)] <- NA
			       	y
			       }))
		if (log && is.null(zero))
			data[[type]][data[[type]] == 0] <- NA
		data$pty <-
		if (type == "interval")
			factor(sign(data$delta - delta),
			       levels = c(0, -1, 1),
			       labels = c("main", "short", "long"))
		else factor("main")

		## Axis limits (x)
		xlim <- xlim_bak
		if (is.null(xlim))
			xlim <- range(data$time)

		## Axis limits (y)
		ylim <- ylim_bak
		if (is.null(ylim)) {
			y <- data[[type]]
			if (type == "rt" && show_predict > 0L)
				y <- c(y, unlist1(cache_predict[elu]))
			y <- y[!is.na(y)]
			if (length(y) == 0L || all(y == 0))
				ylim <- c(if (log) 0.1 else 0, 1)
			else if (!log)
				ylim <- c(0, 1.04 * max(y))
			else if (type == "rt") {
				ry <- range(y)
				ylim <- exp(base::log(ry) +
				            c(-1, 1) * 0.04 * diff(base::log(ry)))
				ylim[!is.finite(ylim)] <- ry[!is.finite(ylim)]
			}
			else if (all(y == 0 | y >= 1) && any(y > 1))
				ylim <- exp(c(-0.04, 1.04) * log(max(y)))
			else {
				ry <- range(y[y > 0])
				ylim <- exp(base::log(ry) +
				            c(-1, 1) * 0.04 * diff(base::log(ry)))
				ylim[!is.finite(ylim)] <- ry[!is.finite(ylim)]
			}
		}
		if (log && !is.null(zero))
			data[[type]][data[[type]] == 0] <-
				if (is.na(zero)) ylim[1L] else zero

		plot.new()
		plot.window(xlim = xlim, ylim = ylim,
		            xaxs = "i", yaxs = "i",
		            log = if (log) "y" else "")

		## Fitting windows
		if (is.list(args <- control$window)) {
			args[c("xleft", "xright")] <- data_windows[c("start", "end")]
			args[c("ybottom", "ytop")] <- as.list(ylim)
			do.call(rect, args)
		}

		## Observed data
		for (s in levels(data$pty))
			if (is.list(args <- control$data[[s]])) {
				reserved <- c("formula", "data", "subset", "x", "y")
				args[names(args) %in% reserved] <- NULL
				args <- c(list(formula = formula, data = data,
				               subset = (data$pty == s)),
				          args)
				do.call(points, args)
			}

		## Confidence bands on predicted curves
		if (show_predict == 2L && is.list(args <- control$predict$ci))
			for (px in split(cache_predict, cache_predict$window, drop = TRUE)) {
				if (type == "cumulative") {
					c0 <- data[[type]][match(px$time[1L], data$time, 0L)]
					if (is.na(c0))
						next
					px$lower <- c0 + px$lower - px$lower[1L]
					px$upper <- c0 + px$upper - px$upper[1L]
				}
				args$x <- c(px$time, rev(px$time))
				args$y <- c(px$lower, rev(px$upper))
				do.call(polygon, args)
			}

		## Predicted curves
		if (show_predict > 0L && is.list(args <- control$predict$value))
			for (px in split(cache_predict, cache_predict$window, drop = TRUE)) {
				if (type == "cumulative") {
					c0 <- data[[type]][match(px$time[1L], data$time, 0L)]
					if (is.na(c0))
						next
					px$value <- c0 + px$value - px$value[1L]
				}
				reserved <- c("formula", "data", "subset", "x", "y")
				args[names(args) %in% reserved] <- NULL
				args <- c(list(formula = value ~ time, data = px), args)
				do.call(lines, args)
			}

		## Asymptote
		if (show_asymptote > 0L && is.list(args <- control$asymptote)) {
			args[c("x0", "x1")] <- data_windows[c("start", "end")]
			args[c("y0", "y1")] <- cache_r["value"]
			do.call(segments, args)
		}

		## Box
		if (is.list(args <- control$box))
			do.call(box, args)

		## Axis (x)
		if (is.list(args <- control$axis$x)) {
			args$side <- 1
			do.call(switch(time_as, Date = Daxis, baxis), args)
		}

		## Axis title (x)
		if (time_as == "numeric" && is.list(args <- control$title$xlab)) {
			args$xlab <- xlab
			do.call(title, args)
		}

		## Axis (y)
		if (is.list(args1 <- control$axis$y)) {
			args1$side <- 2
			args1$las <- 1
			args1$at <- axTicks(side = 2)
			args1$labels <-
				if (type != "rt" && max(args1$at) >= 1e+05)
					get_scientific_labels(args1$at)
				else as.character(args1$at)
			do.call(baxis, args1)
			width <- max(strwidth(args1$labels,
			                      units = "inches",
			                      cex = args1$cex.axis,
			                      font = args1$font.axis,
			                      family = args1$family))
			line <- args1$mgp[2L] +
				diff(grconvertX(c(0, width), "inches", "lines")) + 0.75
		}

		## Axis title (y)
		if (is.list(args2 <- control$title$ylab)) {
			args2$ylab <- ylab[[1L]]
			if (is.list(args1))
				args2$line <- line
			else {
				line <- args2$line
				if (is.null(line))
					line <- args2$line <- args2$mgp[1L]
			}
			do.call(title, args2)
			width <- strheight(args2$ylab,
			                   units = "inches",
			                   cex = args2$cex.lab,
			                   font = args2$font.lab,
			                   family = args2$family)
			line <- line +
				diff(grconvertX(c(0, width), "inches", "lines")) + 1.25
		}

		if (type == "rt" && log) {
			## Axis (y), outer
			if (is.list(args1)) {
				par(new = TRUE)
				plot.window(xlim = xlim, ylim = base::log(2) / ylim,
				            xaxs = "i", yaxs = "i", log = "y")
				args1$at <- axTicks(side = 2)
				args1$labels <- as.character(args1$at)
				args1$mgp <- line + args1$mgp
				do.call(baxis, args1)
				width <- max(strwidth(args1$labels,
				                      units = "inches",
				                      cex = args1$cex.axis,
				                      font = args1$font.axis,
				                      family = args1$family))
				line <- args1$mgp[2L] +
					diff(grconvertX(c(0, width), "inches", "lines")) + 0.75
				par(new = TRUE)
				plot.window(xlim = xlim, ylim = ylim,
				            xaxs = "i", yaxs = "i", log = "y")
			}

			## Axis title (y), outer
			if (is.list(args2)) {
				args2$ylab <- ylab[[2L]]
				args2$line <- line
				do.call(title, args2)
			}
		}

		## Initial doubling times
		if (show_doubling > 0L && is.list(args <- control$doubling)) {
			doubling <- as.list(base::log(2) / cache_r[elu])
			names(doubling) <- elu[c(1L, 3L, 2L)]

			doubling$labels <- c(ci = sprintf("(%.3g%% CI)", 100 * level),
			                     value = "value",
			                     legend = "initial doubling time:")
			w <-
			function(s, l) {
				if (!is.list(l))
					return(NA_real_)
				strwidth(s,
				         units = "user",
				         cex = l$cex * gp$cex,
				         font = l$font,
				         family = l$family)
			}
			doubling$widths <- mapply(w,
			                          s = doubling$labels,
			                          l = args[names(doubling$labels)])
			h <-
			function(s, l) {
				if (!is.list(l))
					return(NA_real_)
				height <- strheight(s,
				                    units = "inches",
				                    cex = l$cex * gp$cex,
				                    font = l$font,
				                    family = l$family)
				diff(grconvertY(c(0, height), "inches", "lines"))
			}
			doubling$heights <- mapply(h,
			                           s = doubling$labels,
			                           l = args[names(doubling$labels)])

			doubling$lines <- rep.int(0.25, 5L)
			names(doubling$lines) <-
				paste0(rep.int(c("", "label_"), c(2L, 3L)),
				       names(doubling$labels)[c(1:2, 1:3)])

			if (show_doubling == 2L && is.list(args$ci))
				doubling$lines[-1L] <- doubling$lines[["ci"]] +
					doubling$heights[["ci"]] + 0.15
			if (is.list(args$value))
				doubling$lines[-(1:2)] <- doubling$lines[["value"]] +
					doubling$heights[["value"]] + 1
			if (show_doubling == 2L && is.list(args$ci))
				doubling$lines[-(1:3)] <- doubling$lines[["label_ci"]] +
					doubling$heights[["ci"]] + 0.15
			if (is.list(args$value))
				doubling$lines[-(1:4)] <- doubling$lines[["label_value"]] +
					doubling$heights[["value"]] + 0.25

			## Legend
			if (show_legend <- is.list(args$legend)) {
				adj <- args$legend$adj
				x_legend <- xlim[1L] + adj *
					max(0, xlim[2L] - xlim[1L] - doubling$widths[["legend"]])
				x_body <- x_legend + doubling$widths[["legend"]] -
					0.5 * max(doubling$widths[c("value", "ci")])

				args$legend$text <- doubling$labels[["legend"]]
				args$legend$side <- 3
				args$legend$line <- doubling$lines[["label_legend"]]
				args$legend$at <- x_legend
				args$legend$adj <- 0
				args$legend$padj <- 0
				args$legend$cex <- args$legend$cex * gp$cex
				do.call(mtext, args$legend)
			}

			## Values
			if (is.list(args$value)) {
				args$value$text <-
					c(sprintf("%.1f", doubling$value),
					  if (show_legend) doubling$labels[["value"]])
				args$value$side <- 3
				args$value$line <-
					c(rep.int(doubling$lines[["value"]],
					          nrow(data_windows)),
					  if (show_legend) doubling$lines[["label_value"]])
				args$value$at <-
					c((data_windows$start + data_windows$end) / 2,
					  if (show_legend) x_body)
				args$value$adj <- 0.5
				args$value$padj <- 0
				args$value$cex <- args$value$cex * gp$cex
				do.call(mtext, args$value)
			}

			## Confidence intervals
			if (show_doubling == 2L && is.list(args$ci)) {
				args$ci$text <-
					c(sprintf("(%.1f, %.1f)", doubling$lower, doubling$upper),
					  if (show_legend) doubling$labels[["ci"]])
				args$ci$side <- 3
				args$ci$line <-
					c(rep.int(doubling$lines[["ci"]], nrow(data_windows)),
					  if (show_legend) doubling$lines[["label_ci"]])
				args$ci$at <-
					c((data_windows$start + data_windows$end) / 2,
					  if (show_legend) x_body)
				args$ci$adj <- 0.5
				args$ci$padj <- 0
				args$ci$cex <- args$ci$cex * gp$cex
				do.call(mtext, args$ci)
			}
		}

		## Plot subtitle
		if (is.list(args1 <- control$title$sub)) {
			args1$main <- sub[i]
			names(args1) <- base::sub("\\.sub$", ".main", names(args1))
			if (show_doubling > 0L && is.list(control$doubling))
				args1$line <- doubling$lines[["value"]] +
					doubling$heights[["value"]] + 1
			else if (is.null(args1$line))
				args1$line <- args1$mgp[1L] + 1
			do.call(title, args1)
			height <- strheight(args1$main,
			                    units = "inches",
			                    cex = args1$cex.main,
			                    font = args1$font.main,
			                    family = args1$family)
			line <- args1$line +
				diff(grconvertY(c(0, height), "inches", "lines")) + 0.25
		}

		## Plot title
		if (is.list(args2 <- control$title$main)) {
			args2$main <- main[i]
			if (is.list(args1))
				args2$line <- line
			else if (show_doubling > 0L && is.list(control$doubling))
				args2$line <- doubling$lines[["value"]] +
					doubling$heights[["value"]] + 1
			do.call(title, args2)
		}
	} # loop over plots

	invisible(NULL)
}
