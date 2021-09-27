## Trying to be defensive... in the unlikely situation that the namespace
## is not already loaded, then loading it with 'loadNamepace' inside of
## a test has the side effect of registering S3 methods. This side effect
## is not reversed by 'unloadNamespace', so it would persist during later
## tests...
skip_if_not("epigrowthfit" %in% loadedNamespaces())

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

test_that("version mismatch", {
  path_to_build_versions <- system.file("build_versions", package = "epigrowthfit")
  expect_true(nzchar(path_to_build_versions))

  build_versions <- readLines(path_to_build_versions)
  expect_gt(length(build_versions), 0L)
  expect_true(all(grepl("^[[:alnum:].]{2,} [[:digit:]]+(\\.[[:digit:]]+)+$", build_versions)))
  expect_true(all(vapply(strsplit(build_versions, " "), `[[`, "", 1L) %in% rownames(installed.packages())))

  v <- packageVersion("TMB")
  mmp <- unclass(v)[[1L]]
  vm <- as.package_version("0.0.0")
  vp <- as.package_version(as.numeric_version(list(c(mmp[1:2], mmp[3L] + 1L))))

  test_versions <- paste("TMB", c(v, vm, vp))
  regexp <- list(NA, NULL, NULL)
  unloadNamespace("epigrowthfit")
  for (i in 1:3) {
    writeLines(test_versions[[i]], path_to_build_versions)
    eval(bquote(expect_warning(loadNamespace("epigrowthfit"), regexp = .(regexp[[i]]))))
    unloadNamespace("epigrowthfit")
  }
  writeLines(build_versions, path_to_build_versions)
  loadNamespace("epigrowthfit")
})
