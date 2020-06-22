## setwd("../data")
ff <- list.files(pattern="\\.*[Rr][dD].*")
xoresave <- function(f) {
    L <- load(f)
    save(list=L,file=f,version=2)
}
lapply(ff,resave)
