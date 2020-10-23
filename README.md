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
version 3.5.0 or greater and R package
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
compatible with their version of R. Rtools is needed to
build from source R packages that use C++ code, and
**epigrowthfit** is such a package. Once installed, Rtools
must be added to the search path. Users with Rtools 40
(compatible with R Version 4.0.0 or greater) can accomplish
this by running

```r
write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"',
  file = "~/.Renviron",
  append = TRUE,
  sep = "\n"
)
```

then restarting R, as explained
[here](https://cran.r-project.org/bin/windows/Rtools/).
Older versions of Rtools are typically installed to
`C:\Rtools`. Users with one of these versions can run

```r
write('PATH="C:\\Rtools\\bin;${PATH}"',
  file = "~/.Renviron",
  append = TRUE,
  sep = "\n"
)
```

replacing the substring `C:\\Rtools` if necessary.
Finally, at least one Windows user running Rtools 40
has encountered the following install-time error:

```
C:/Rtools/mingw_64/bin/g++: No such file or directory
```

Compiler `g++` was not found because `C:\Rtools` was not
the path to Rtools 40. The error was resolved after running

```r
write('BINPATH="${RTOOLS40_HOME}\\mingw${WIN}\\bin"',
  file = "~/.Renviron",
  append = TRUE,
  sep = "\n"
)
```

and restarting R, as this specified the correct path
to the compiler.

### Installing from the remote repository

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

### Installing from a local repository

If you have cloned **epigrowthfit** from GitHub, then you can install
it on the command line from your local repository. Local installation
requires R package
[**devtools**](https://CRAN.R-project.org/package=devtools):

```r
if (!require(devtools)) {
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
library(epigrowthfit)

## This vignette
vignette("epigrowthfit-vignette")

## Help pages
?"epigrowthfit-package" # package
?data_set_name          # data set "data_set_name"
?function_name          # function "function_name"
?"class_name-methods"   # S3 methods for class "class_name"
```
