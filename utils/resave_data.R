## setwd("../data")
ff <- list.files(pattern = "\\.[Rr][dD]a(ta)?$")
resave <- function(f) {
    L <- load(f)
    save(list = L, file = f, version = 2, compress = "xz")
}
lapply(ff, resave)
