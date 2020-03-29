## install all packages needed by epigrowthfit
pkgs <- c("bbmle","emdbook","Hmisc", "Deriv", "lubridate", "zoo", "numDeriv",
          "Rdpack","testthat", "outbreaks","knitr")
i1 <- installed.packages()
pkgs <- pkgs[!pkgs %in% rownames(i1)]
install.packages(pkgs)




