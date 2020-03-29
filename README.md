
# epigrowthfit

`epigrowthfit` is an R package for estimating parameters associated with initial epidemic growth.  The methodology is based on:

[Ma J, Dushoff J, Bolker BM, Earn DJD (2014). “Estimating initial epidemic growth rates.” Bull. Math. Biol., 76(1), 245-260.](https://davidearn.mcmaster.ca/publications/MaEtAl2014)

## installation

- run `remotes::install_github("davidearn/epigrowthfit")` to install the `epigrowthfit` package
- this should automatically install required upstream packages, but if necessary you can download and run [install_pkgs.R](./install_pkgs.R) (please [add a comment here](https://github.com/davidearn/epigrowthfit/issues/1) if you have to run this step, or for other installation problems)

## notes

- We have chosen to make the current package available publicly now (March 2020) in order to help people who are analyzing COVID-19 data.  The documentation will improve over time.  If you have suggestions for improvement, please e-mail us (`bolker@mcmaster.ca` and `earn@math.mcmaster.ca`).
- The current defaults and settings were developed for analysis of historical plague epidemics; in particular, any estimates of R0 and final size **are based on generation intervals for bubonic plague**; we are working to disable output of this information for now


