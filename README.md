# epigrowthfit

**epigrowthfit** is an R package for estimating parameters associated with
epidemic growth. Methods are adapted from Ma *et al*. (2014) and Earn *et al*
(2020):

[Ma J, Dushoff J, Bolker BM, Earn DJD (2014). “Estimating initial epidemic growth rates.” *Bulletin of Mathematical Biology*, **76**, 245-260.](https://davidearn.mcmaster.ca/publications/MaEtAl2014)

Earn DJD, Ma J, Poinar HN, Dushoff J, Bolker BM (2020). “Acceleration of plague outbreaks in the second pandemic.” *Proceedings of the National Academy of Sciences USA*, in press.

## Package installation

### Dependencies

**epigrowthfit** depends on an installation of
[R](https://www.r-project.org/)
version 3.3.0 or greater and R package
[**TMB**](https://CRAN.R-project.org/package=TMB),
and imports from R packages
[**mathjaxr**](https://CRAN.R-project.org/package=mathjaxr),
[**Rdpack**](https://CRAN.R-project.org/package=Rdpack),
and
[**emdbook**](https://CRAN.R-project.org/package=emdbook).
Building the vignette requires additional R packages
[**knitr**](https://CRAN.R-project.org/package=knitr)
and
[**shape**](https://CRAN.R-project.org/package=shape),
as well as a [LaTeX](https://www.latex-project.org/)
distribution.

Windows users require an installation of
[Rtools](https://cran.r-project.org/bin/windows/Rtools/)
compatible with their version of R. Once installed, Rtools
must be added to the search path. With Rtools 40, this can
be done by running

```r
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
```

then restarting R, as explained
[here](https://cran.r-project.org/bin/windows/Rtools/).
With older versions of Rtools, the following line of code
should work instead, assuming that Rtools was installed to
`C:\Rtools`.

```r
writeLines('PATH="C:\\Rtools\\bin;${PATH}"', con = "~/.Renviron")
```

### From the remote repository

**epigrowthfit** can be installed from this GitHub repository
using function `install_github()` from the
[**remotes**](https://CRAN.R-project.org/package=remotes)
package:

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

To disable automatic installation of dependencies,
set `dependencies = FALSE`. To avoid building the
vignette, which requires a LaTeX distribution and
can add a few minutes to the installation time if
many LaTeX packages must also be installed, set
`dependencies = NA` and `build_vignettes = FALSE`.

### From a local repository

If you have cloned **epigrowthfit** from GitHub, then you can install
it on the command line from your local repository. Locall installation
requires R package
[**devtools**](https://CRAN.R-project.org/package=devtools):

```r
if (!require(devtools)) {
  install.packages("devtools")
}
```

Rnter the `devel` branch of the repository with
`git checkout devel`, then run `make` in the root directory,
namely `epigrowthfit/`. This will install **epigrowthfit**
and its dependencies and build the vignette. You can verify
that the installation was successful by running
`(require(epigrowthfit))` in R.

## Documentation

Package documentation can be accessed as follows:

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
