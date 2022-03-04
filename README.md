# epigrowthfit

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/davidearn/epigrowthfit/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/davidearn/epigrowthfit/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**epigrowthfit** is an R package for fitting nonlinear mixed effects
models of epidemic growth to collections of one or more disease
incidence time series. It can be applied to birth processes other
than epidemics, as the statistical machinery is agnostic to the
precise interpretation of supplied count data. **epigrowthfit**
is built on [Template Model Builder](https://cran.r-project.org/package=TMB).

## Installation

**epigrowthfit** is not (yet) available on
[CRAN](https://cran.r-project.org/),
but can be installed from its sources on GitHub.

```r
remotes::install_github("davidearn/epigrowthfit", build_vignettes = TRUE)
```

Since the package contains compiled code, installation from source 
depends on compilers and related tools.
These will already be available on most modern Linux installations.
Windows users must have installed
[Rtools](https://cran.r-project.org/bin/windows/Rtools/)
and placed it on their `PATH`.
macOS users must have installed Apple's Command Line Tools for 
[Xcode](https://developer.apple.com/xcode/).
The most recent version of Command Line Tools can be installed 
by running 

```shell
sudo rm -rf /Library/Developer/CommandLineTools
sudo xcode-select --install
```

in Terminal. Older binaries can be downloaded
[here](https://developer.apple.com/download/all/?q=xcode).

Vignette builds depend on a 
[LaTeX](https://www.latex-project.org/get/) distribution.
Specifically, `PATH` must specify a directory containing `pdflatex`. 
A minimal distribution can be installed and managed via the R package 
[**tinytex**](https://cran.r-project.org/package=tinytex).
Should errors due to vignette builds persist, one can always install
**epigrowthfit** without vignettes by setting `build_vignettes = FALSE`.

Installation-related issues should be reported
[here](https://github.com/davidearn/epigrowthfit/issues/1).

### OpenMP support

**epigrowthfit** supports low level parallelism via 
[OpenMP](https://en.wikipedia.org/wiki/OpenMP).
Whether OpenMP is supported on a given system depends on:

1. the platform;
2. the compiler and compiler flags used to build R from sources
   (these settings are preserved in `R_HOME/etc/Makeconf`, 
   unless it has been edited manually); and
3. the compiler and compiler flags used to build **epigrowthfit**
   from sources (these settings are taken from `HOME/.R/Makevars`
   [`HOME/.R/Makevars.win` on Windows] and `R_HOME/etc/Makeconf`, 
   in that order).

On modern Linux and Windows (with Rtools), if one installs
a current version of R from a binary prebuilt by CRAN, then
OpenMP support is often automatic.
On macOS, some configuration is needed, because Apple `clang`,
the C compiler provided with Apple's Command Line Tools,
does not support OpenMP.
Below are _provisional_ instructions for obtaining a `clang`
toolchain supporting OpenMP on Macs with Intel- or ARM-based 
architectures running native builds of R. 
These may need to be adapted to future releases of macOS, R,
and R prerequisites.

First, install Command Line Tools:

```shell
sudo xcode-select --install
```

Next, find your version of Apple `clang` with

```shell
clang --version
```

Use [this](https://mac.r-project.org/openmp/) link to download 
the OpenMP runtime library suitable for your macOS and `clang` 
versions, and install by unpacking to root.
For example, if your `clang` version is `clang-1205.0.22.11`, 
then you could do:

```shell
wget https://mac.r-project.org/openmp/openmp-11.0.1-darwin20-Release.tar.gz
sudo tar xvf openmp-11.0.1-darwin20-Release.tar.gz -C /
```

Finally, add these lines to `HOME/.R/Makevars`:

```make
CPPFLAGS+=-I/usr/local/include -Xclang -fopenmp
LDFLAGS+=-L/usr/local/lib -lomp
```

### Package version mismatch

**epigrowthfit** must be binary-compatible with 
[**TMB**](https://cran.r-project.org/package=TMB). 
In turn, **TMB** must be binary-compatible with 
[**Matrix**](https://cran.r-project.org/package=Matrix). 
This means that the version of **TMB** currently installed should 
match the version that was installed when **epigrowthfit** was 
built from sources. 
Similarly, the version of **Matrix** currently installed should match
the version that was installed when **TMB** was built from sources. 
Thus, version mismatch can occur when a user:

* updates **TMB** without rebuilding **epigrowthfit** from source, or
* updates **Matrix** without rebuilding **TMB** from source.

**epigrowthfit** and **TMB** perform checks when loaded 
(e.g., via `library`) and issue warnings if they detect a mismatch. 
If you are warned about a **Matrix**-**TMB** mismatch, then run:

```r
install.packages("TMB", type = "source")
remotes::install_github("davidearn/epigrowthfit")
```

If you are only warned about a **TMB**-**epigrowthfit** mismatch,
then only run the second line.

## References

[Kristensen K, Nielsen A, Berg CW, Skaug H, Bell BM (2016). "TMB: Automatic differentiation and Laplace approximation." *Journal of Statistical Software*, **70**, 1-21.](https://www.jstatsoft.org/article/view/v070i05)

[Ma J, Dushoff J, Bolker BM, Earn DJD (2014). “Estimating initial epidemic growth rates.” *Bulletin of Mathematical Biology*, **76**, 245-260.](https://davidearn.mcmaster.ca/publications/MaEtAl2014)

[Earn DJD, Ma J, Poinar HN, Dushoff J, Bolker BM (2020). “Acceleration of plague outbreaks in the second pandemic.” *Proceedings of the National Academy of Sciences USA*, **117**, 27703-27711.](https://davidearn.mcmaster.ca/publications/EarnEtAl2020)
