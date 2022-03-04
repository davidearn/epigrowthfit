test_that("basic", {
    o <- egf_cache("egf-1.rds")
    capture.output({
        expect_condition(print(o), regexp = NA)
        expect_invisible(print(o))
        expect_identical(print(o), o)
    })
})
