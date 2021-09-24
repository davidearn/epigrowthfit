# epigrowthfit

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/davidearn/epigrowthfit/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/davidearn/epigrowthfit/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/davidearn/epigrowthfit/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/davidearn/epigrowthfit/actions/workflows/test-coverage.yaml)
<!-- badges: end -->

**epigrowthfit** is an R package for fitting nonlinear mixed effects
models of epidemic growth to collections of one or more disease
incidence time series. It can be applied to birth processes other
than epidemics, as the statistical machinery is agnostic to the
precise interpretation of supplied count data. **epigrowthfit**
is built on [Template Model Builder](https://github.com/kaskr/adcomp).

## Installation

**epigrowthfit** is not (yet) available on
[CRAN](https://cran.r-project.org/). 
It can be installed from its sources on GitHub.

```r
remotes::install_github("davidearn/epigrowthfit")
```

Since the package contains compiled code, installation from source 
depends on compilers and related tools.
These will already be available on most modern Linux installations.
Windows users must have installed
[Rtools](https://cran.r-project.org/bin/windows/Rtools/)
and added it to their `PATH`.
macOS users must have installed Apple's Command Line Tools,
which are packaged with [Xcode](https://developer.apple.com/xcode/).
The most recent version of Xcode can be installed by running 
`xcode-select --install` in a terminal.
However, users should install whichever version of Xcode was used 
to build their R binary.
General information about installing R packages, 
including extensive platform-specific notes, can be found in `R-admin`
(i.e., the R manual called _R Installation and Administration_,
which can be browsed using `help.start()`).

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
On macOS, some configuration is needed, because the default
Command Line Tools `clang` toolchain does not include OpenMP
support. _Provisional_ instructions for obtaining a working
`clang` toolchain on Macs with Intel- or ARM-based architectures,
running _native_ builds of R, are as follows.

First, install Xcode.

```bash
xcode-select --install
```

Next, install the LLVM `clang` toolchain, which includes OpenMP 
support:

```bash
## With Homebrew
brew update
brew install llvm

## Without Homebrew (Intel-based Macs only)
wget https://github.com/llvm/llvm-project/releases/download/llvmorg-11.0.0/clang+llvm-11.0.0-x86_64-apple-darwin.tar.xz
tar xvf clang+llvm-11.0.0-x86_64-apple-darwin.tar.xz -C /usr/local --strip-components 1
rm clang+llvm-11.0.0-x86_64-apple-darwin.tar.xz
```

Finally, add these lines to `HOME/.R/Makevars` (creating the 
file if necessary and deleting those lines not indicated for 
your system):

```make
#### Intel-based Macs only ####

R_HOME=/usr/local
LLVM_HOME=/usr/local/opt/llvm

## Mojave and later
SDK_PATH=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
CC="${LLVM_HOME}/bin/clang -isysroot ${SDK_PATH}"
CXX="${LLVM_HOME}/bin/clang++ -isysroot ${SDK_PATH}"

## High Sierra and earlier
CC=${LLVM_HOME}/bin/clang
CXX=${LLVM_HOME}/bin/clang++

## Xcode 12 and later
CFLAGS="-g -O2 -Wall -pedantic -Wno-implicit-function-declaration"
CXXFLAGS="-g -O2 -Wall -pedantic"

## Xcode 11 and earlier
CFLAGS="-g -O2 -Wall -pedantic"
CXXFLAGS="-g -O2 -Wall -pedantic"


#### ARM-based Macs only ####

R_HOME=/opt/R/arm64
LLVM_HOME=/opt/homebrew/opt/llvm
SDK_PATH=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk

CC="${LLVM_HOME}/bin/clang -isysroot ${SDK_PATH} -target arm64-apple-macos11"
CXX="${LLVM_HOME}/bin/clang++ -isysroot ${SDK_PATH} -target arm64-apple-macos11"

CFLAGS="-falign-functions=8 -g -O2 -Wall -pedantic -Wno-implicit-function-declaration"
CXXFLAGS="-g -O2 -Wall -pedantic"


#### Both Intel- and ARM-based Macs ####

SHLIB_OPENMP_CFLAGS=-fopenmp
SHLIB_OPENMP_CXXFLAGS=-fopenmp

CPPFLAGS="-I${LLVM_HOME}/include -I${R_HOME}/include"
LDFLAGS="-L${LLVM_HOME}/lib -L${R_HOME}/lib"
```

These steps, gathered from different sections of `R-admin`, 
should ensure that building **epigrowthfit** from source 
results in an installation that supports OpenMP parallelism. 
If this is not the case, then please comment on 
[#1](https://github.com/davidearn/epigrowthfit/issues/1)
with your system details so that these instructions can be 
updated.


### Package version mismatch

**epigrowthfit** must be binary-compatible with **TMB**. 
In turn, **TMB** must be binary-compatible with **Matrix**. 
This means that the version of **TMB** currently installed should 
match the version that was installed when **epigrowthfit** was 
built from source. Similarly, the version of **Matrix** currently 
installed should match the version that was installed when 
**TMB** was built from source. Thus, version mismatch can occur
when a user:

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
then the second line should suffice.

## References

[Kristensen K, Nielsen A, Berg CW, Skaug H, Bell BM (2016). "TMB: Automatic differentiation and Laplace approximation." *Journal of Statistical Software*, **70**, 1-21.](https://www.jstatsoft.org/article/view/v070i05)

[Ma J, Dushoff J, Bolker BM, Earn DJD (2014). “Estimating initial epidemic growth rates.” *Bulletin of Mathematical Biology*, **76**, 245-260.](https://davidearn.mcmaster.ca/publications/MaEtAl2014)

[Earn DJD, Ma J, Poinar HN, Dushoff J, Bolker BM (2020). “Acceleration of plague outbreaks in the second pandemic.” *Proceedings of the National Academy of Sciences USA*, **117**, 27703-27711.](https://davidearn.mcmaster.ca/publications/EarnEtAl2020)
