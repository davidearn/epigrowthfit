.onLoad <- function(libname, pkgname)
	library.dynam("epigrowthfit", pkgname, libname)

.onUnload <- function(libpath)
	library.dynam.unload("epigrowthfit", libpath)
