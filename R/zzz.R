#' @importFrom utils packageVersion
.onLoad <- function(libname, pkgname) {
    build_versions <- system.file("build_versions", package = pkgname,
                                  mustWork = TRUE)
    for (l in strsplit(readLines(build_versions), " ")) {
        depname <- l[1L]
        buildver <- l[2L]
        currentver <- as.character(packageVersion(depname))
        if (buildver != currentver) {
            warning(wrap("Package version mismatch detected. ",
                         paste(sQuote(depname), currentver), " is installed, ",
                         "but ", sQuote(pkgname), " was built with ",
                         paste(sQuote(depname), buildver), ". ",
                         "Install ", paste(sQuote(depname), buildver), " or ",
                         "a version of ", sQuote(pkgname), " built with ",
                         paste(sQuote(depname), currentver), "."))
        }
    }
    ## FIXME: Can this line replace 'useDynLib' directive in NAMESPACE?
    ## library.dynam("epigrowthfit", pkgname, libname)
    invisible(NULL)
}

.onUnload <- function(libpath) {
    library.dynam.unload("epigrowthfit", libpath)
    invisible(NULL)
}
