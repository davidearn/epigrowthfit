#' @importFrom utils packageVersion
.onLoad <- function(libname, pkgname) {
  l <- strsplit(readLines(system.file("build_versions", package = "epigrowthfit")), " ")
  for (i in seq_along(l)) {
    depname <- l[[i]][[1L]]
    build <- l[[i]][[2L]]
    current <- as.character(packageVersion(depname))
    if (build != current) {
      warning(wrap(
        "Package version mismatch detected. ",
        paste(sQuote(depname), current), " is installed, ",
        "but ", sQuote(pkgname), " was built with ", paste(sQuote(depname), build), ". ",
        "Install ", paste(sQuote(depname), build), " or ",
        "a version of ", sQuote(pkgname), " built with ", paste(sQuote(depname), current), "."
      ))
    }
  }
  invisible(NULL)
}
