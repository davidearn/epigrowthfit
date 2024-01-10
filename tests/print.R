library(epigrowthfit)
options(warn = 2L, error = recover)


## basic ######
    o <- egf_cache("egf-1.rds")
    capture.output({
        expect_condition(print(o), regexp = NA)
        expect_invisible(print(o))
        identical(print(o), o)
    })

