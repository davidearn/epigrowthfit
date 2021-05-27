gradwrap <- function(gr, inner.control=NULL, inner.method=NULL) {
  function(...) {
    g1 <- gr(...)
    if (!(length(g1)==1 && is.na(g1))) return(g1) ## OK?
    ee <- environment(gr)
    orig_ic <- ee$inner.control
    orig_im <- ee$inner.method
    on.exit( {
      ee$inner.control <- orig_ic
      ee$inner.method  <- orig_im
    })
    ## etc.
    g2 <- gr(...)
    if (<...>) stop("too bad, I really tried")
    return(g2)
  }
