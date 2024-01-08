print.egf <- function(x, width = 0.9 * getOption("width"), indent = 2L, ...) {
    width <- as.integer(width)
    indent <- strrep(" ", indent)


    ## Set up ==========================================================

    ## Top level nonlinear model
    mu <- quote(f(t) - f(s))
    if (x$model$day_of_week) {
        mu <- call("*", quote(w(s, t)), mu)
    }
    e1 <- call("~",
               quote(X(s, t)),
               switch(x$model$family,
                      pois = call("dpois", lambda = mu),
                      nbinom = call("dnbinom", mu = mu, size = quote(disp))))
    e2 <- call("=",
               quote(f(t)),
               switch(x$model$curve,
                      exponential = quote(c0 * exp(r * t)),
                      subexponential = quote(c0 * (1 + alpha * (1 - p) * t / c0^(1 - p))^(1 / (1 - p))),
                      gompertz = quote(K * exp(-exp(-alpha * (t - tinfl)))),
                      logistic = quote(K / (1 + exp(-r * (t - tinfl)))),
                      richards = quote(K / (1 + a * exp(-a * r * (t - tinfl)))^(1 / a))))
    if (x$model$excess) {
        e2[[3L]] <- call("+", quote(b * t), e2[[3L]])
    }
    lines_top <- c(deparse(e1), "where", deparse(e2))

    ## Bottom level mixed effects model
    line <- function(s, offset) {
        formula <- call("~", str2lang(s), formula(x, top = s)[[2L]])
        paste0(strrep(" ", offset), paste0(deparse(formula), collapse = "\n"))
    }
    offset <- function(s) max(n <- nchar(s)) - n
    names_top <- egf_get_names_top(x, link = TRUE)
    lines_bottom <- mapply(line, names_top, offset(names_top))

    ## Number of observations
    frame <- model.frame(x)
    n <- c(sum(!is.na(frame$x)), nlevels(frame$window), nlevels(frame$ts))
    c1 <- sprintf("%d", n)
    c2 <- c("observation", "fitting window", "time series")
    c2[1:2] <- pluralize(c2[1:2], n[1:2])
    lines_nobs <- align(c1, c2, justify = c("r", "l"))


    ## Write to 'stdout'  ==============================================

    heading("Top level nonlinear model", width = width)
    cat("\n")
    writeLines(paste0(indent, lines_top))
    cat("\n")
    str(x$model, no.list = TRUE, indent.str = indent, give.head = FALSE)
    cat("\n")

    heading("Bottom level mixed effects model", width = width)
    cat("\n")
    writeLines(paste0(indent, lines_bottom))
    cat("\n")

    heading("Data", width = width)
    cat("\n")
    writeLines(paste0(indent, lines_nobs))

    invisible(x)
}

print.egf_no_fit <- print.egf
