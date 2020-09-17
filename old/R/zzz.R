.onLoad <- function(libname, pkgname) {
    ## prevent NOTE false positives from R CMD check ...
    N <- a <- b <- k <- x <- y <- NULL

    ## export derivative rules (not sure this is necessary ...)

    drule[["lbeta"]] <- drule[["w_lbeta"]] <- alist(a=dfun(a,b),
                                                    b=dfun(b,a))
    drule[["dfun"]] <- alist(x=dfun2(x,y),
                             y=dfun2(y,x))
    
    drule[["NBconst"]] <- alist(x=dfun(x,k), k=dfun(k,x) + 1/x)
}
