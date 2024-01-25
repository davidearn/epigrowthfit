.onLoad <-
function(libname, pkgname) {
	ns <- parent.env(environment())

	.dll <- paste0("epigrowthfit", .Platform[["dynlib.ext"]])
	.dll <-
		if (nzchar(arch <- .Platform[["r_arch"]]))
			system.file("libs", arch, .dll, package = pkgname, mustWork = TRUE)
		else
			system.file("libs",       .dll, package = pkgname, mustWork = TRUE)
	assign(".dll", envir = ns, inherits = FALSE, .dll)

	library.dynam("epigrowthfit", pkgname, libname)
}

.onUnload <-
function(libpath)
	library.dynam.unload("epigrowthfit", libpath)
