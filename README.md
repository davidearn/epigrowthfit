# epigrowthfit

**epigrowthfit** is an R package for fitting nonlinear mixed effects
models of epidemic growth to collections of one or more disease
incidence time series.  It can be applied to birth processes other than
epidemics, as the statistical machinery is agnostic to the precise
interpretation of supplied count data.  **epigrowthfit** is built on
[Template Model Builder](https://cran.r-project.org/package=TMB).

## Installation

[CRAN](https://cran.r-project.org/package=epigrowthfit) distributes
both the package sources and binaries for Windows and macOS.  Hence
typical users will install **epigrowthfit** with

```r
install.packages("epigrowthfit")
```

or perhaps

```r
install.packages("epigrowthfit", type = "source")
```

to force installation from sources where installation of a binary
would occur by default.  The rest of this section concerns the
`type = "source"` case.

Installation from sources depends on compilers and related tools.
These will already be available on modern Linux installations.
Windows users must have installed
[Rtools](https://cran.r-project.org/bin/windows/Rtools/).
macOS users must have installed Apple's Command Line Tools for
[Xcode](https://developer.apple.com/xcode/)
and GNU Fortran.
The most recent version of Command Line Tools supporting the
version of macOS in use can be installed by running

```shell
sudo rm -rf /Library/Developer/CommandLineTools
sudo xcode-select --install
```

in Terminal.
Binaries for older versions of Command Line Tools can be downloaded
[here](https://developer.apple.com/download/all/?q=Command%20Line%20Tools%20for%20Xcode).
GNU Fortran should be installed following the instructions on the
R for macOS Developers [page](https://mac.r-project.org/tools/).

Vignette building depends on a
[LaTeX](https://www.latex-project.org/get/) distribution.
Specifically, `PATH` must specify the location of `pdflatex`.

Issues related to installation should be reported (with relevant output)
[here](https://github.com/davidearn/epigrowthfit/issues/1).
Notably, details for Windows and macOS can differ for non-current
versions of R or non-standard installations of R, where
the "standard" way to install R is to download and unpack a binary
built and published by [CRAN](https://cran.r-project.org/).

Compiler errors encountered on macOS are almost always explained
by unmet dependencies or masking of native tools, headers, and
libraries with non-native ones (e.g., ones installed by Homebrew).
Masking occurs due to dubious configuration of `PATH` or dubious
setting (typically in `~/.R/Makevars`) of Make variables such as
`CPPFLAGS` and `LDFLAGS`.  Users should reattempt compilation
after removing suspicious components of `PATH` (e.g., by
removing relevant lines of startup files in your home directory,
then launching a new shell) and (re)moving `~/.R/Makevars`.

**epigrowthfit** supports low level parallelism via
[OpenMP](https://en.wikipedia.org/wiki/OpenMP).  The support tends to
be automatic except on macOS, because Apple's Command Line Tools do
not bundle the requisite headers and runtime library.  The R for macOS
Developers [page](https://mac.r-project.org/openmp/) makes these files
available for download, hence macOS users can follow the installation
instructions there.  If `clang --version` gives 1403.0.22.14.1 (say),
then one would do

```shell
curl -O https://mac.r-project.org/openmp/openmp-15.0.7-darwin20-Release.tar.gz
sudo mkdir -p /opt/R/$(uname -m)
sudo tar -xvf openmp-15.0.7-darwin20-Release.tar.gz --strip-components=2 -C /opt/R/$(uname -m)
```

and create a `~/.R/Makevars` containing the lines

```make
CPPFLAGS += -Xclang -fopenmp
LDFLAGS += -lomp
```

The flags in the `tar` command line ensure that the files are unpacked
under `/opt/R/x86_64` (Intel) or `/opt/R/arm64` (Apple Silicon).
Standard installations of R will already be configured to search there
for dependencies.

## Repository structure

Active development happens on branch `master`.  Tested changes intended
for the next release are ported to branch `release-candidate`,
where tarballs submitted to CRAN are eventually built.  Neither `master`
nor `release-candidate` should be considered stable.

The stable branches are named `release-x.y.z`.  They branch from
`release-candidate` before the version number there is incremented,
typically just after a tarball is submitted to CRAN.

To install **epigrowthfit** from sources in a given branch or commit,
install [**remotes**](https://cran.r-project.org/package=remotes) and
run, e.g.,

```r
remotes::install_github("davidearn/epigrowthfit", ref = "release-0.15.2")
remotes::install_github("davidearn/epigrowthfit", ref = "cf6fdd8")
```

## References

[Kristensen K, Nielsen A, Berg CW, Skaug H, Bell BM (2016). "TMB: Automatic differentiation and Laplace approximation." *Journal of Statistical Software*, **70**, 1-21.](https://www.jstatsoft.org/article/view/v070i05)

[Ma J, Dushoff J, Bolker BM, Earn DJD (2014). “Estimating initial epidemic growth rates.” *Bulletin of Mathematical Biology*, **76**, 245-260.](https://davidearn.mcmaster.ca/publications/MaEtAl2014)

[Earn DJD, Ma J, Poinar HN, Dushoff J, Bolker BM (2020). “Acceleration of plague outbreaks in the second pandemic.” *Proceedings of the National Academy of Sciences USA*, **117**, 27703-27711.](https://davidearn.mcmaster.ca/publications/EarnEtAl2020)
