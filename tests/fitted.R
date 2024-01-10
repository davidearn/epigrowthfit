library(epigrowthfit)
options(warn = 2L, error = recover)


## object ######
    o <- egf_cache("egf-1.rds")
    fo <- egf_cache("fitted-egf-1.rds")
    fo_expected <-
        data.frame(top = gl(2L, 20L, labels = c("log(r)", "log(c0)")),
                   ts = gl(10L, 2L, 40L, labels = LETTERS[1:10]),
                   window = gl(20L, 1L, 40L,
                               labels = sprintf("window_%02d", 1:20)),
                   estimate = as.double(o$sdreport$value),
                   se = as.double(o$sdreport$sd))
    attr(fo_expected, "se") <- TRUE
    class(fo_expected) <- c("egf_fitted", "data.frame")
    identical(fo, fo_expected)


## confint ######
    fo <- egf_cache("fitted-egf-1.rds")
    cfo <- confint(fo, level = 0.95)
    cfo_expected <- fo
    cfo_expected[c("lower", "upper")] <- wald(fo$estimate, fo$se, level = 0.95)
    cfo_expected[["se"]] <- NULL
    attr(cfo_expected, "level") <- 0.95
    attr(cfo_expected, "se") <- NULL
    class(cfo_expected) <- "data.frame"
    identical(cfo, cfo_expected)

