ff <- list.files()
resave <- function(f) {
    L <- load(f)
    save(list=L,file=f,version=2)
}
lapply(ff,resave)
