library(epigrowthfit)
options(warn = 2L, error = recover)


## egf_get_names_top ######
    x <- egf_get_names_top(NULL, link = FALSE)
    expect_type(x, "character")
    expect_gt(length(x), 0L)
    expect_false(anyNA(x))
    expect_named(x, NULL)

    model <- egf_model()
    x0 <- egf_get_names_top(model, link = FALSE)
    expect_type(x0, "character")
    expect_gt(length(x0), 0L)
    expect_true(all(x0 %in% x))

    o <- list(model = model)
    class(o) <- "egf"
    identical(egf_get_names_top(o, link = FALSE), x0)


## egf_has_random ######
    e <- new.env()
    e$data <- list()
    o <- list(tmb_out = list(env = e))
    class(o) <- "egf"

    e$data$Z <- matrix(double(9L), 3L, 3L)
    expect_true(egf_has_random(o))

    e$data$Z <- matrix(double(0L), 3L, 0L)
    expect_false(egf_has_random(o))


## egf_has_converged ######
    o <- list(optimizer_out = list(convergence = 0L),
              value = 100,
              gradient = runif(10, -0.5, 0.5),
              hessian = TRUE)
    class(o) <- c("egf", "list")

    expect_true(egf_has_converged(o))
    expect_false(egf_has_converged(within(o, optimizer_out$convergence <- 1L)))
    expect_false(egf_has_converged(within(o, value <- NaN)))
    expect_false(egf_has_converged(within(o, gradient[1L] <- 10)))
    expect_false(egf_has_converged(within(o, hessian <- FALSE)))
    identical(egf_has_converged(within(o, hessian <- NA)), NA)


## egf_(expand|condense)_par ######
    o <- egf_cache("egf-2.rds")
    par <- o$tmb_out$env$last.par.best
    len <- c(table(factor(names(par), levels = c("beta", "theta", "b"))))
    epar <- egf_expand_par(o$tmb_out, par = par)
    all.equal(epar, structure(c(par[1:4], theta = 0, par[-(1:4)]),
                                 lengths = len + c(0L, 1L, 0L)))
    cpar <- egf_condense_par(o$tmb_out, par = epar)
    all.equal(cpar, structure(par, lengths = len))


## egf_get_sdreport ######
    o <- egf_cache("egf-1.rds")
    sdr <- o$sdreport
    identical(egf_get_sdreport(o)$cov.fixed, sdr$cov.fixed)
    o$sdreport <- NULL
    expect_warning(all.equal(egf_get_sdreport(o)$cov.fixed, sdr$cov.fixed))


## egf_preprofile ######
    ## FIXME: Only testing a trivial case as no elements of 'beta' are mapped
    o <- egf_cache("egf-1.rds")
    l <- egf_preprofile(o, subset = seq_len(20L), top = c("log(r)", "log(c0)"))
    expect_type(l, "list")
    expect_length(l, 2L)
    expect_named(l, c("Y", "A"), ignore.order = TRUE)
    identical(unname(l$Y), Matrix(0, 20L, 2L))
    identical(unname(l$A), sparseMatrix(i = seq_len(40L),
                                               j = rep.int(1:2, c(20L, 20L)),
                                               x = 1,
                                               dims = c(40L, 5L)))


## egf_cache ######
    subdir <- tools::R_user_dir("epigrowthfit", "cache")
    lf1 <- list.files(subdir)
    file <- "test.rds"
    a <- 1
    x <- egf_cache(file, a)
    y <- egf_cache(file, a <- 2)
    lf2 <- list.files(subdir)
    z <- egf_cache(file, clear = TRUE)
    lf3 <- list.files(subdir)

    identical(x, 1)
    identical(y, 1)
    identical(z, 0L)
    identical(a, 1)
    identical(setdiff(lf2, lf1), file)
    identical(lf3, lf1)

