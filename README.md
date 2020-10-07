# epigrowthfit

**epigrowthfit** is an R package for estimating parameters associated with
initial epidemic growth. Methods are adapted from Ma *et al*. (2014):

[Ma J, Dushoff J, Bolker BM, Earn DJD (2014). “Estimating initial epidemic growth rates.” Bull. Math. Biol., 76(1), 245-260.](https://davidearn.mcmaster.ca/publications/MaEtAl2014)

## Installation

### Local

To install from a local repository:

* Check that you are in the `devel` branch with `git branch`.
* Check that you are in the root directory of the package.
* Run `make`.

### Remote

To install from the remote repository, run this script in R:

```r
if (!require(remotes)) install.packages("remotes")
remotes::install_github("davidearn/fastbeta", ref = "devel", build_vignettes = TRUE)
```

## Documentation

Package documentation can be browsed like so:

```r
library(epigrowthfit)

## List of exported functions and data sets
ls("package:epigrowthfit")

## Vignette
vignette("epigrowthfit-vignette")

## Help pages
#options(help_type = "html") # if running R from command line
?"epigrowthfit-package"   # package
?function_name            # exported function "function_name"
?"class_name-methods"     # S3 methods for class "class_name"
```
