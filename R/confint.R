confint.egf <- function(object,
                        parm,
                        level = 0.95,
                        top = egf_get_names_top(object, link = TRUE),
                        link = TRUE,
                        method = c("wald", "profile", "uniroot"),
                        parallel = egf_parallel(),
                        trace = FALSE,
                        grid_len = 12,
                        interval_scale = 7,
                        breaks = NULL,
                        probs = NULL,
                        subset = NULL,
                        append = NULL,
                        .subset = NULL,
                        .append = NULL,
                        ...) {
    stopifnot(is_number_in_interval(level, 0, 1, "()"), is_true_or_false(link))
    method <- match.arg(method)
    elu <- c("estimate", "lower", "upper")

    ## Initial doubling time and basic reproduction number are
    ## monotone functions of initial exponential growth rate,
    ## so confidence intervals can be obtained "for free" _if_
    ## the top level nonlinear model _has_ an initial exponential
    ## growth rate
    names_top <- names_top_aug <- egf_get_names_top(object, link = TRUE)
    mono <- c("tdoubling", "R0")
    if (object$model$curve %in% c("exponential", "logistic", "richards")) {
        names_top_aug <- c(names_top_aug, mono)
    }
    top <- top_aug <- unique(match.arg(top, names_top_aug, several.ok = TRUE))
    if ("R0" %in% top) {
        stopifnot(!is.null(breaks), !is.null(probs))
    }
    top[top %in% mono] <- "log(r)"
    top <- unique(top)

    frame_windows <- model.frame(object, "windows")
    frame_combined <- model.frame(object, "combined")
    subset <- if (is.null(.subset)) substitute(subset) else .subset
    subset <- egf_eval_subset(subset, frame_combined, parent.frame())
    append <- if (is.null(.append)) substitute(append) else .append
    append <- egf_eval_append(append, frame_combined, baseenv())

    if (method == "wald") {
        fo <- fitted(object, top = top, se = TRUE,
                     .subset = subset, .append = append)
        res <- confint(fo, level = level, link = link)

    } else if (method == "profile") {
        po <- profile(object,
                      level = level + min(0.01, 0.1 * (1 - level)),
                      top = top,
                      parallel = parallel,
                      trace = trace,
                      grid_len = grid_len,
                      .subset = subset,
                      .append = append)
        res <- confint(po, level = level, link = link)
        res[["linear_combination"]] <- NULL
        attr(res, "A") <- attr(res, "x") <- NULL

    } else { # "uniroot"
        stopifnot(inherits(parallel, "egf_parallel"),
                  is_true_or_false(trace),
                  is_number(interval_scale, "positive"))
        n <- sum(!object$random)

        l <- egf_preprofile(object, subset = subset, top = top)
        Y <- l$Y
        A <- l$A

        m <- nrow(A)
        a <- lapply(seq_len(m), function(i) A[i, ])

        ## y := nll_restricted - nll_minimum = 0.5 * deviance
        target <- 0.5 * qchisq(level, df = 1)
        sd.range <- interval_scale
        nomp <- object$control$omp_num_threads

        do_uniroot <- function(i, a) {
            if (trace) {
                cat(sprintf("Computing confidence interval %d of %d...\n",
                            i, m))
            }
            res <- TMB::tmbroot(obj, lincomb = a, target = target,
                                sd.range = sd.range, trace = FALSE)
            names(res) <- c("lower", "upper")
            res
        }

        if (parallel$method == "snow") {
            environment(do_uniroot) <- .GlobalEnv

            ## Reconstruct list of arguments to 'MakeADFun'
            ## from object internals for retaping
            args <- egf_tmb_remake_args(object$tmb_out, par = object$best)

            ## Retrieve path to shared object for loading
            dll <- system.file("libs", TMB::dynlib("epigrowthfit"),
                               package = "epigrowthfit", mustWork = TRUE)

            cl <- parallel$cl
            if (is.null(cl)) {
                cl <- do.call(makePSOCKcluster, parallel$args)
                on.exit(stopCluster(cl), add = TRUE)
            }
            vars <- c("dll", "nomp", "args", "trace", "m", "target", "sd.range")
            clusterExport(cl, varlist = vars, envir = environment())
            clusterEvalQ(cl, {
                dyn.load(dll)
                if (TMB::openmp(n = NULL) > 0L) {
                    TMB::openmp(n = nomp)
                }
                obj <- do.call(TMB::MakeADFun, args)
            })
            res <- clusterMap(cl, do_uniroot, i = seq_len(m), a = a)
        } else {
            if (given_outfile <- nzchar(parallel$outfile)) {
                outfile <- file(parallel$outfile, open = "wt")
                sink(outfile, type = "output")
                sink(outfile, type = "message")
                on.exit(add = TRUE, {
                    sink(type = "message")
                    sink(type = "output")
                })
            }
            onomp <- TMB::openmp(n = NULL)
            if (onomp > 0L) {
                TMB::openmp(n = nomp)
                on.exit(TMB::openmp(n = onomp), add = TRUE)
            }
            obj <- object$tmb_out
            res <- switch(parallel$method,
                          multicore = do.call(mcMap, c(list(f = do_uniroot, i = seq_len(m), a = a), parallel$args)),
                          serial = Map(do_uniroot, i = seq_len(m), a = a))
        }

        res <- data.frame(top = rep(factor(top, levels = names_top),
                                    each = length(subset)),
                          ts = frame_windows$ts[subset],
                          window = frame_windows$window[subset],
                          estimate = as.double(A %*% object$best[!object$random]),
                          do.call(rbind, res),
                          frame_combined[subset, append, drop = FALSE],
                          row.names = NULL,
                          check.names = FALSE,
                          stringsAsFactors = FALSE)
        res[elu] <- res[elu] + as.double(Y)

        if (!link) {
            f <- lapply(egf_link_extract(levels(res$top)), egf_link_match,
                        inverse = TRUE)
            res[elu] <- in_place_ragged_apply(res[elu], res$top, f = f)
            levels(res$top) <- egf_link_remove(levels(res$top))
        }
    } # "uniroot"

    if (any(mono %in% top_aug)) {
        s <- if (link) "log(r)" else "r"
        res_r <- res[res$top == s, , drop = FALSE]
        if (link) {
            res_r[elu] <- exp(res_r[elu])
        }

        if ("tdoubling" %in% top_aug) {
            eul <- elu[c(1L, 3L, 2L)]
            res_tdoubling <- res_r
            res_tdoubling[elu] <- log(2) / res_r[eul]
            res_tdoubling$top <- factor("tdoubling")
            res <- rbind(res, res_tdoubling)
        }

        if ("R0" %in% top_aug) {
            res_R0 <- res_r
            res_R0[elu] <- lapply(res_r[elu], compute_R0,
                                  breaks = breaks, probs = probs)
            res_R0$top <- factor("R0")
            res <- rbind(res, res_R0)
        }

        if (!"log(r)" %in% top_aug) {
            res <- res[res$top != s, , drop = FALSE]
        }
    }

    row.names(res) <- NULL
    attr(res, "method") <- method
    attr(res, "level") <- level
    attr(res, "frame_windows") <- frame_windows # for 'plot.egf_confint'
    class(res) <- c("egf_confint", oldClass(res))
    res
}

plot.egf_confint <- function(x,
                             type = c("bars", "boxes"),
                             time_as = c("Date", "numeric"),
                             per_plot = switch(type, bars = 12L, boxes = 3L),
                             subset = NULL,
                             order = NULL,
                             main = NULL,
                             label = NULL,
                             ...) {
    type <- match.arg(type)
    time_as <- match.arg(time_as)
    stopifnot(is_number(per_plot, "positive", integer = TRUE))
    per_plot <- as.integer(per_plot)

    subset <- egf_eval_subset(substitute(subset), x, parent.frame())
    if (length(subset) == 0L) {
        stop1("'subset' indexes zero confidence intervals, ",
              "so there is nothing to plot.")
    }
    label <- egf_eval_label(substitute(label), x, parent.frame())
    order <- egf_eval_order(substitute(order), x, parent.frame())
    subset <- order[order %in% subset]

    a <- attributes(x)
    if (is.null(main)) {
        s <- switch(type, bars = "fitting window", boxes = "time series")
        main <- sprintf("%.3g%% confidence intervals by %s", 100 * a$level, s)
    }
    if (is.null(label)) {
        s <- switch(type, bars = "window", boxes = "ts")
        label <- as.character(x[[s]])
    }

    nx <- c("top", "ts", "window", "estimate", "lower", "upper")
    x <- x[nx]
    x$label <- label
    x <- x[subset, , drop = FALSE]
    x$top <- factor(x$top)

    if (type == "bars") {
        plot.egf_confint.bars(x, per_plot = per_plot, main = main)
    } else {
        i <- match(x$window, a$frame_windows$window, 0L)
        frame_windows <- a$frame_windows[i, c("start", "end"), drop = FALSE]
        x <- data.frame(x, frame_windows)
        plot.egf_confint.boxes(x, per_plot = per_plot, main = main,
                               time_as = time_as)
    }
}

plot.egf_confint.bars <- function(x, per_plot, main) {
    x <- split(x, x$top)
    for (xlab in names(x)) { # loop over parameters
        data <- x[[xlab]]
        n <- nrow(data)
        argna <- lapply(data[c("lower", "upper")], is.na)
        xlim <- range(data[c("estimate", "lower", "upper")], na.rm = TRUE)

        i <- 0L
        while (i < n) { # loop over plots
            k <- i + seq_len(min(per_plot, n - i))
            plot.new()
            plot.window(xlim = xlim, ylim = c(per_plot + 1, 0),
                        xaxs = "r", yaxs = "i")
            gp <- par(c("usr", "cex", "mar"))
            sfcex <- get_sfcex(x$label,
                               target = 0.95 * max(0, gp$mar[2L] - 0.25),
                               units = "lines")
            abline(v = axTicks(side = 1), lty = 3, col = "grey75")
            segments(x0 = replace(data$lower[k], argna$lower[k], gp$usr[1L]),
                     x1 = data$estimate[k],
                     y0 = seq_along(k),
                     y1 = seq_along(k),
                     lty = 1 + as.double(argna$lower[k]),
                     lwd = 2 - as.double(argna$lower[k]))
            segments(x0 = data$estimate[k],
                     x1 = replace(data$upper[k], argna$upper[k], gp$usr[2L]),
                     y0 = seq_along(k),
                     y1 = seq_along(k),
                     lty = 1 + as.double(argna$upper[k]),
                     lwd = 2 - as.double(argna$upper[k]))
            points(x = data$estimate[k],
                   y = seq_along(k),
                   pch = 21,
                   bg = "grey80")
            box()
            axis(side = 1)
            mtext(text = data$label[k],
                  side = 2,
                  line = 0.25,
                  at = seq_along(k),
                  las = 1,
                  adj = 1,
                  padj = 0.5,
                  cex = gp$cex * sfcex)
            title(xlab = xlab)
            title(main, adj = 0)
            i <- i + per_plot
        } # loop over plots
    } # loop over parameters

    invisible(NULL)
}

plot.egf_confint.boxes <- function(x, per_plot, main, time_as) {
    xlim <- range(x[c("start", "end")])
    x <- split(x, x$top)

    op <- par(mfrow = c(per_plot, 1))
    on.exit(par(op))

    for (ylab in names(x)) { # loop over parameters
        data <- x[[ylab]]
        ylim <- range(data[c("estimate", "lower", "upper")], na.rm = TRUE)
        data <- split(data, factor(data$ts, levels = unique(data$ts)))
        n <- length(data)

        i <- 0L
        while (i < n) { # loop over pages
            for (k in i + seq_len(min(per_plot, n - i))) { # loop over plots
                plot.new()
                plot.window(xlim = xlim, ylim = ylim)
                gp <- par(c("mfrow", "usr", "mai", "omi", "mgp",
                            "cex", "cex.axis", "cex.lab"))

                v <- Daxis(side = 1, major = NULL, minor = NULL)$minor
                abline(v = v, lty = 3, col = "grey75")

                argna <- lapply(data[[k]][c("lower", "upper")], is.na)
                data[[k]]$lower[argna$lower] <- gp$usr[3L]
                data[[k]]$upper[argna$upper] <- gp$usr[4L]

                for (l in seq_len(nrow(data[[k]]))) { # loop over CI
                    rect(xleft = data[[k]]$start[l],
                         xright = data[[k]]$end[l],
                         ybottom = data[[k]]$lower[l],
                         ytop = data[[k]]$estimate[l],
                         col = if (argna$lower[l]) NA else "grey80",
                         border = "grey50",
                         lty = if (argna$lower[l]) 2 else 1)
                    rect(xleft = data[[k]]$start[l],
                         xright = data[[k]]$end[l],
                         ybottom = data[[k]]$estimate[l],
                         ytop = data[[k]]$upper[l],
                         col = if (argna$upper[l]) NA else "grey80",
                         border = "grey50",
                         lty = if (argna$upper[l]) 2 else 1)
                    segments(x0 = data[[k]]$start[l],
                             x1 = data[[k]]$end[l],
                             y0 = data[[k]]$estimate[l],
                             y1 = data[[k]]$estimate[l],
                             col = "grey50",
                             lwd = 2)
                } # loop over CI

                p <- diff(grconvertY(c(0, 0.08), "npc", "inches"))
                px <- diff(grconvertX(c(0, p), "inches", "user"))
                py <- diff(grconvertY(c(0, p), "inches", "user"))

                par(cex = 1)
                text(x = gp$usr[1L] + px,
                     y = gp$usr[4L] - py,
                     labels = data[[k]]$label[1L],
                     adj = c(0, 1),
                     cex = gp$cex.lab)
                par(cex = gp$cex)

                box()
                axis(side = 2)
            } # loop over plots

            par(cex = 1)
            if (time_as == "numeric") {
                axis(side = 1)
            } else {
                Daxis(side = 1,
                      major = list(cex.axis = 1.2 * gp$cex.axis,
                                   mgp = gp$mgp + c(0, 1.5, 0),
                                   tick = FALSE))
            }
            par(cex = gp$cex)

            ## Hack to avoid dealing with normalized figure coordinates

            par(new = TRUE, mfrow = c(1, 1),
                mai = gp$mai + gp$omi, omi = c(0, 0, 0, 0))
            plot.window(xlim = xlim, ylim = c(0, 1))
            title(main = main, adj = 0)
            title(ylab = ylab, adj = 0.5)
            par(new = FALSE, mfrow = gp$mfrow,
                mai = gp$mai, omi = gp$omi)

            i <- i + per_plot
        } # loop over pages
    } # loop over parameters

    invisible(NULL)
}
