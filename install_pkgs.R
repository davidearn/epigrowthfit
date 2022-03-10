## Install all 'epigrowthfit' dependencies
pkgs <- c("Deriv", "numDeriv", "bbmle", "emdbook", "Hmisc", "zoo", "lubridate",
          "Rdpack", "knitr", "tikzDevice", "dplyr", "tidyr", "ggplot2",
          "outbreaks")
i1 <- installed.packages()
install.packages(setdiff(pkgs, rownames(i1)))
