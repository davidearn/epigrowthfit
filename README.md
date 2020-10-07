# epigrowthfit

This branch houses **epigrowthfitPNAS**, an alias for the version of R
package **epigrowthfit** used to support analysis of plague epidemics
in Earn *et al*. (2020):

The sole purpose of this branch is to enable simultaneous installation
of this version and the development version of **epigrowthfit** housed
in branch `devel`, as the development version is not fully backwards
compatible.

## Installation

To install, run:

```r
remotes::install_github("davidearn/epigrowthfit",
                        ref = "pnas",
			dependencies = TRUE,
			build_vignettes = TRUE)			
```


