# epigrowthfit

**epigrowthfit** is an R package for fitting nonlinear mixed
effects models of epidemic growth to disease incidence time
series. Methods are adapted from Ma *et al*. (2014) and Earn
*et al* (2020):

[Ma J, Dushoff J, Bolker BM, Earn DJD (2014). “Estimating initial epidemic growth rates.” *Bulletin of Mathematical Biology*, **76**, 245-260.](https://davidearn.mcmaster.ca/publications/MaEtAl2014)

[Earn DJD, Ma J, Poinar HN, Dushoff J, Bolker BM (2020). “Acceleration of plague outbreaks in the second pandemic.” *Proceedings of the National Academy of Sciences USA*, **117**(44), 27703-27711.](https://davidearn.mcmaster.ca/publications/EarnEtAl2020)

## Package installation

### Dependencies

**epigrowthfit** depends on an installation of
[R](https://www.r-project.org/)
3.5.0 or greater and R package
[**TMB**](https://CRAN.R-project.org/package=TMB),
and imports from R packages
[**Matrix**](https://CRAN.R-project.org/package=Matrix)
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
compatible with their version of R. Rtools is needed to
build from source R packages like **epigrowthfit** that
rely on compiled C++ code. Once installed, Rtools must be
added to `PATH`. Users with Rtools 40
(compatible with R 4.0.0 or greater) can accomplish
this by running

```r
write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"',
  file = "~/.Renviron",
  append = TRUE,
  sep = "\n"
)
```

then restarting R. Users with Rtools 35 (compatible with
R 3.3.0&ndash;3.6.3) must set both `PATH` and `BINPREF` by running

```r
write('PATH="C:\\Rtools\\bin;${PATH}"',
  file = "~/.Renviron",
  append = TRUE,
  sep = "\n"
)
write('BINPREF="C:\\Rtools\\mingw_${WIN}\\bin"',
  file = "~/.Renviron",
  append = TRUE,
  sep = "\n"
)
```

then restarting R. The above assumes that Rtools 35 was installed to
`C:\Rtools` (typically the default location).


### Installing from GitHub

**epigrowthfit** can be installed from this GitHub repository
using function `install_github()` from the
[**remotes**](https://CRAN.R-project.org/package=remotes)
package:

```r
if (!require("remotes")) {
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

### Installing from a local repository

If you have cloned **epigrowthfit** from GitHub, then you can install
it on the command line from your local repository. Local installation
requires R package
[**devtools**](https://CRAN.R-project.org/package=devtools):

```r
if (!require("devtools")) {
  install.packages("devtools")
}
```

Enter the `devel` branch of the repository with
`git checkout devel`, then run `make` in the root directory,
namely `epigrowthfit/`. This will install **epigrowthfit**
and its dependencies and build the vignette. You can verify
that the installation was successful by running
`(require(epigrowthfit))` in R.

## Documentation

Package documentation can be accessed as follows:

```r
library("epigrowthfit")

## Vignette
vignette("epigrowthfit-vignette")

## Help pages
?"epigrowthfit-package"   # package
?data_set_name            # data set "data_set_name"
?function_name            # function "function_name"
?generic_name.class_name  # S3 methods for class "class_name"
```
