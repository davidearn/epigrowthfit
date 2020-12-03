get_wave <- function(n, subset) {
  subset <- check_subset(n, subset)
  wave <- integer(n)
  for (i in seq_along(subset)) {
    wave[subset[[i]]] <- i
  }
  factor(wave, levels = seq_along(subset))
}


