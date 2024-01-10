## skip_if_not(is.null(pkgload::dev_meta("epigrowthfit")))

test_that(".on(Load|Unload)", {
    is_ns_loaded <- function() "epigrowthfit" %in% loadedNamespaces()
    is_so_loaded <- function() "epigrowthfit" %in% names(getLoadedDLLs())
    unloadNamespace("epigrowthfit")
    expect_false(is_ns_loaded())
    expect_false(is_so_loaded())
    loadNamespace("epigrowthfit")
    expect_true(is_ns_loaded())
    expect_true(is_so_loaded())
})
