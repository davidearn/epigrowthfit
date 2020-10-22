# epigrowthfit

**epigrowthfit** is an R package for estimating parameters associated with
initial epidemic growth. Methods are adapted from Ma *et al*. (2014)
and Earn *et al* (2020):

[Ma J, Dushoff J, Bolker BM, Earn DJD (2014). “Estimating initial epidemic growth rates.” *Bulletin of Mathematical Biology*, **76**, 245-260.](https://davidearn.mcmaster.ca/publications/MaEtAl2014)

Earn DJD, Ma J, Poinar HN, Dushoff J, Bolker BM (2020). “Acceleration of plague outbreaks in the second pandemic.” *Proceedings of the National Academy of Sciences USA*, in press.

## Installation

### Dependencies

**epigrowthfit** depends on installation of
[R](https://www.r-project.org/) version 3.3.0 or greater and
R package [**TMB**](https://CRAN.R-project.org/package=TMB).
**epigrowthfit** imports from R packages
[**mathjaxr**](https://CRAN.R-project.org/package=mathjaxr)
(for LaTeX in help pages),
[**Rdpack**](https://CRAN.R-project.org/package=Rdpack)
(for bibliographic references in help pages), and
[**emdbook**](https://CRAN.R-project.org/package=emdbook)
(for evaluation of the Lambert W function).
Compilation of the **epigrowthfit** vignette requires
additional R packages
[**knitr**](https://CRAN.R-project.org/package=knitr) and
[**shape**](https://CRAN.R-project.org/package=shape),
as well as a full distribution of
[LaTeX](https://www.latex-project.org/).
Windows users require an installation of
[Rtools](https://cran.r-project.org/bin/windows/Rtools/)
compatible with their version of R and must add Rtools
to their search path. For R version 4.0.0 or greater,
this can be done by running

```r
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
```

then restarting R, as explained
[here](https://cran.r-project.org/bin/windows/Rtools/).
For older R versions, this can be done by running

```r
writeLines('PATH="C:\\Rtools\\bin;${PATH}"', con = "~/.Renviron")
```

then restarting R, assuming that Rtools was installed to `C:\Rtools`.

### Remote

To install from the remote repository, run this script in R:

```r
if (!require(remotes)) {
  install.packages("remotes")
}
remotes::install_github("davidearn/epigrowthfit",
  ref = "devel",
  dependencies = TRUE,
  build_vignettes = TRUE
)
```

### Local

To install from a local repository:

* Check that you are in the `devel` branch with `git branch`.
* Check that you are in the root directory of the package.
* Run `make`.

## Documentation

Package documentation can be browsed like so:

```r
library(epigrowthfit)

## This vignette
vignette("epigrowthfit-vignette")

## Help pages
?"epigrowthfit-package" # package
?data_set_name          # data set "data_set_name"
?function_name          # function "function_name"
?"class_name-methods"   # S3 methods for class "class_name"
```
