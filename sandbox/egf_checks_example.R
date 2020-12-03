source("R/utils.R")
source("R/egf_utils.R")
source("R/egf_checks.R")

date <- rep(as.character(as.Date(0:99, origin="1995-01-01")), 4)
cases <- sample(0:1e03, size = 400, replace = TRUE)
continent <- factor(unlist(lapply(c("europe", "asia"), rep, 200)))
country <- factor(unlist(lapply(c("france", "germany", "china", "japan"), rep, 100)))
index <- integer(400) # fitting windows
index[40:70] <- 1L
index[140:170] <- 2L
index[240:270] <- 3L
index[340:370] <- 4L
index <- factor(index, exclude = 0L)
data <- data.frame(date, cases, continent, country, index)

formula <- cases ~ date
fixed <- list(r = ~country, p = ~continent)
random <- NULL
par_names <- c("r", "thalf", "K", "p")
dfmt <- "%Y-%m-%d"

formula <- check_formula(formula)
fixed <- check_fixed(fixed, par_names)
random <- check_random(random, par_names)

frame <- check_data(formula, data, index, fixed, random, dfmt)
madf_data <- make_madf_data(frame, FALSE, FALSE)


