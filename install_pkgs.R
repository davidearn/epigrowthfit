## install all packages needed by epigrowthfit
pkgs <- c("bbmle","emdbook","Hmisc", "Deriv", "lubridate", "zoo", "numDeriv",
          "Rdpack","testthat", "outbreaks","knitr")
pkgs <- pkgs[!pkgs %in% rownames(i1)]
install.packages(pkgs)




