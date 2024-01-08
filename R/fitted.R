fitted.egf <- function(object,
                       top = egf_get_names_top(object, link = TRUE),
                       link = TRUE,
                       se = FALSE,
                       subset = NULL,
                       append = NULL,
                       .subset = NULL,
                       .append = NULL,
                       ...) {
    stopifnot(is_true_or_false(link),
              is_true_or_false(se))
    if (se && !link) {
        stop1("Standard errors are not available for inverse link-transformed ",
              "fitted values.")
    }

    names_top <- egf_get_names_top(object, link = TRUE)
    top <- unique(match.arg(top, names_top, several.ok = TRUE))

    frame_windows <- model.frame(object, "windows")
    frame_combined <- model.frame(object, "combined")
    subset <- if (is.null(.subset)) substitute(subset) else .subset
    subset <- egf_eval_subset(subset, frame_combined, parent.frame())
    append <- if (is.null(.append)) substitute(append) else .append
    append <- egf_eval_append(append, frame_combined, baseenv())

    if (se) {
        sdr <- egf_get_sdreport(object)
        ssdr <- summary(sdr, select = "report")
        index <- rownames(ssdr) == "Y"
        Y <- ssdr[index, "Estimate"]
        Y_se <- ssdr[index, "Std. Error"]
        dim(Y) <- dim(Y_se) <- object$tmb_out$env$ADreportDims$Y
    } else {
        Y <- object$tmb_out$report(object$best)$Y
    }

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
                      frame_combined[subset, append, drop = FALSE],
                      row.names = NULL,
                      check.names = FALSE)
    attr(res, "se") <- se
    class(res) <- c("egf_fitted", oldClass(res))
    res
}

fitted.egf_no_fit <- function(object,
                              top = egf_get_names_top(object, link = TRUE),
                              link = TRUE,
                              se = FALSE,
                              subset = NULL,
                              append = NULL,
                              .subset = NULL,
                              .append = NULL,
                              ...) {
    if (se) {
        stop1("Standard errors cannot be computed until the model ",
              "is estimated. Retry after doing, e.g., ",
              "'object <- update(object, se = TRUE, fit = TRUE, ...)'.")
    }

    ## Passing arguments to method for class "egf" without evaluating
    ## 'subset' or 'append' requires minor acrobatics
    call <- match.call(expand.dots = FALSE)
    call[[1L]] <- quote(fitted.egf)
    call$... <- NULL
    nms <- names(call)
    i <- match(nms[-1L], c("subset", "append"), 0L) == 0L
    call[-1L][i] <- lapply(nms[-1L][i], as.name)

    object$best <- object$init
    eval(call)
}

confint.egf_fitted <- function(object, parm, level = 0.95, link = TRUE, ...) {
    if (!isTRUE(attr(object, "se"))) {
        stop1("'object' must supply link scale fitted values ",
              "and corresponding standard errors. Retry with ",
              "'object = fitted(<\"egf\" object>, link = TRUE, se = TRUE)'.")
    }
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
