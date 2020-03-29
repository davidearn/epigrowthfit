## install all packages needed by epigrowthfit
pkgs <- c("broom.mixed", "colorspace", "cowplot", "Deriv", "dfoptim", "directlabels", "dplyr", "emdbook",
	"emmeans", "epigrowthfit", "ggplot2", "ggstance", "glmmTMB", "gtable",
	"Hmisc", "huxtable", "lubridate", "multcomp", "nloptr", "numDeriv", "optimx", "plyr", "purrr", "RColorBrewer",
        "Rdpack", "tibble", "tidyr", "tikzDevice", "xtable", "remotes", "viridis", "zoo")
i1 <- installed.packages()
pkgs <- pkgs[!pkgs %in% rownames(i1)]
install.packages(pkgs)
library(remotes)
install_github("bbolker/bbmle")




