factor_by_outbreak <- function(data,
                               outbreak_definitions = local({load("outbreak_definitions.RData"); outbreak_definitions})) {
  ## Factor times by outbreak
  data$outbreak <- cut(data$time,
    breaks = c(t(outbreak_definitions[, c("start", "end")])),
    labels = c(rbind(NA, outbreak_definitions$outbreak))[-1],
    right  = FALSE,
    ordered_result = TRUE
  )
  data$outbreak <- droplevels(data$outbreak, exclude = NA)
  data <- data[!is.na(data$outbreak), ]
  ## Factor outbreaks by severity
  data$severity <- data$outbreak
  levels(data$severity) <- outbreak_definitions$severity[match(levels(data$outbreak), outbreak_definitions$outbreak)]
  data$severity <- factor(data$severity, levels = c("minor", "major"), ordered = TRUE)
  data
}
