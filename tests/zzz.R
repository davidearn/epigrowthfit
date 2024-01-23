## .on(Load|Unload) ####################################################

nsLoaded <- function() any("epigrowthfit" == loadedNamespaces())
soLoaded <- function() any("epigrowthfit" == names(getLoadedDLLs()))

loadNamespace("epigrowthfit")
stopifnot(nsLoaded(), soLoaded())
unloadNamespace("epigrowthfit")
stopifnot(!nsLoaded(), !soLoaded())
