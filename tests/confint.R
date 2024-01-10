library(epigrowthfit)
options(warn = 2L, error = recover)


## object ######
    co1 <- egf_cache("confint-egf-1.rds")
    co2 <- egf_cache("confint-egf-2.rds")
    co3 <- egf_cache("confint-egf-3.rds")
    fo <- egf_cache("fitted-egf-1.rds")
    cfo <- confint(fo)
    po <- egf_cache("profile-egf-1.rds")
    cpo <- confint(po)
    cpo[["linear_combination"]] <- NULL

    is.list(co1)
    identical(oldClass(co1), c("egf_confint", "data.frame"))
    all.equal(co1, cfo, ignore_attr = c("class", "method", "frame_windows"))

    is.list(co2)
    identical(oldClass(co2), c("egf_confint", "data.frame"))
    all.equal(co2, cpo,
                 ignore_attr = c("class", "method", "frame_windows", "A", "x"))

    is.list(co3)
    identical(oldClass(co3), c("egf_confint", "data.frame"))
    all.equal(co3, co2, tolerance = 1e-3, ignore_attr = "method")


## parallel ######
    skip_on_cran()
    o <- egf_cache("egf-1.rds")
    co3 <- egf_cache("confint-egf-3.rds")

    co3_multicore <-
        confint(o,
                method = "uniroot",
                subset = (country == "A" & wave == 1),
                parallel = egf_parallel(method = "multicore", cores = 2L))
    all.equal(co3_multicore, co3)

    ## skip_if_not(is.null(pkgload::dev_meta("epigrowthfit")))
    co3_snow <-
        confint(o,
                method = "uniroot",
                subset = (country == "A" & wave == 1),
                parallel = egf_parallel(method = "snow", cores = 2L))
    all.equal(co3_snow, co3)


## plot ######
    skip_on_cran()
    co1 <- egf_cache("confint-egf-1.rds")

    bars <- function() {
        op <- par(mar = c(4.5, 4, 2, 1), oma = c(0, 0, 0, 0))
        on.exit(par(op))
        plot(co1, type = "bars")
    }
    ## vdiffr::expect_doppelganger("plot-egf_confint-1", fig = bars)

    boxes <- function() {
        op <- par(mar = c(0.2, 0, 0.2, 0), oma = c(4.5, 6, 2, 1), las = 1)
        on.exit(par(op))
        plot(co1, type = "boxes")
    }
    ## vdiffr::expect_doppelganger("plot-egf_confint-2", fig = boxes)

    expect_true(TRUE) # otherwise test is considered skipped

