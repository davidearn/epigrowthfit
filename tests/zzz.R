library(epigrowthfit)
options(warn = 2L, error = recover)


## skip_if_not(is.null(pkgload::dev_meta("epigrowthfit")))

## .on(Load|Unload) ######
is_ns_loaded <- function() "epigrowthfit" %in% loadedNamespaces()
is_so_loaded <- function() "epigrowthfit" %in% names(getLoadedDLLs())
unloadNamespace("epigrowthfit")
!is_ns_loaded()
!is_so_loaded()
loadNamespace("epigrowthfit")
is_ns_loaded()
is_so_loaded()

