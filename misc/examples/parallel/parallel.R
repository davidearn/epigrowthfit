library("TMB")
packageVersion("TMB") # 1.7.20
stopifnot(openmp() > 0L)
rds_out <- "res.rds"

rds <- setdiff(list.files(pattern = "\\.rds$"), rds_out)
stopifnot(length(rds) == 1L)
l <- readRDS(rds) # list(data, parameters)

cpp <- list.files(pattern = "\\.cpp$")
stopifnot(length(cpp) > 0L)
for (s in cpp) {
  compile(s, openmp = TRUE)
}

dll <- sub("\\.cpp$", "", cpp)
for (s in dll) {
  dyn.load(dynlib(s))
}

f <- function(dll, n, cores) {
  obj <- MakeADFun(
    data = l$data,
    parameters = l$parameters,
    random = "b",
    DLL = dll,
    inner.method = "newton",
    inner.control = list(maxit = 1000L, trace = 1L),
    silent = FALSE
  )
  benchmark(obj = obj, n = n, cores = cores)
}

res <- sapply(dll, f, n = 1000L, cores = seq_len(min(8L, openmp())), simplify = FALSE)
saveRDS(res, file = rds_out)
