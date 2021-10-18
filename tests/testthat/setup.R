## Example output is cached in package subdirectory 'exdata'.
## Here, we make sure that 'exdata' is populated by running
## all examples whose output we want to test. If examples
## have been run once already (e.g., during 'R CMD check'),
## then running them a second time is harmless, as examples
## also reuse cached output.
topics <- c(
  "egf",
  "confint.egf",
  "fitted.egf",
  "predict.egf",
  "profile.egf",
  "simulate.egf",
  "simulate.egf_model",
  "egf.egf_model_simulate",
  "summary.egf"
)
lapply(topics, example, package = "epigrowthfit", character.only = TRUE, local = TRUE, echo = FALSE)
