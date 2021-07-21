library("TMB")
packageVersion("TMB") # 1.7.20
stopifnot(openmp(NULL) > 0L)

dll <- paste0("parallel", 1:2)
for (s in dll) {
  compile(paste0(s, ".cpp"), openmp = TRUE)
  dyn.load(dynlib(s))
}

a <- 0
b <- 1
log_sigma <- 0

set.seed(151354L)
x <- seq.int(0, 10, length.out = 100001L)
y <- a + b * x + rnorm(x, 0, exp(log_sigma))

data <- list(y = y, x = x)
parameters <- list(a = 0, b = 0, log_sigma = 0)

f <- function(s) {
  MakeADFun(data = data, parameters = parameters, DLL = s)
}
obj <- sapply(dll, f, simplify = FALSE)
res <- lapply(obj, benchmark, n = 1000L, cores = seq_len(min(8L, openmp(NULL))))
saveRDS(res, file = "res.rds")
