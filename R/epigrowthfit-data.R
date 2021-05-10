#' Data in the \pkg{epigrowthfit} package
#'
#' Epidemic time series to which growth rates can be fit,
#' and a few other data sets.
#'
#' @details
#' Below is a list of available data sets with links to their
#' documentation:
#' \describe{
#' \item{[`canadacovid`][canadacovid]}{
#'   Daily confirmations of COVID-19 in Canadian provinces
#'   and territories, from first confirmation to May 8, 2021.
#' }
#' \item{[`covid_generation_interval`][covid_generation_interval]}{
#'   Gamma distribution of the COVID-19 generation interval
#'   fit to data from a cluster of 45 cases in Tianjin, China.
#' }
#' \item{[`husting`][husting]}{
#'   Counts of wills probated in the Court of Husting
#'   during four plague epidemics in 14th century London.
#' }
#' \item{[`canterbury`][canterbury]}{
#'   Counts of wills probated in the Prerogative Court of Canterbury
#'   during 24 plague epidemics in 16th and 17th century London.
#' }
#' \item{[`londonparishes`][londonparishes]}{
#'   Weekly counts of burials listed in extant parish registers
#'   during 24 plague epidemics in 16th and 17th century London.
#' }
#' \item{[`londonbills`][londonbills]}{
#'   Weekly counts of plague deaths recorded in the London Bills of
#'   Mortality during 24 plague epidemics in 16th and 17th century
#'   London.
#' }
#' \item{[`plague_latent_period`][plague_latent_period]}{
#'   Empirical distribution of the latent period of pneumonic plague.
#' }
#' \item{[`plague_infectious_period`][plague_infectious_period]}{
#'   Empirical distribution of the infectious period of pneumonic
#'   plague.
#' }
#' }
#'
#' @name epigrowthfit-data
#' @keywords internal
NULL

#' COVID-19 in Canadian provinces and territories
#'
#' Daily confirmations of COVID-19 in Canadian provinces
#' and territories, from first confirmation to May 8, 2021.
#'
#' @format
#' A data frame with 5451 rows and 4 variables:
#' \describe{
#' \item{`province`}{
#'   A factor. Canadian province or territory (postal abbreviation).
#' }
#' \item{`Date`}{
#'   A Date vector. Daily within each level of `province`,
#'   except prior to invasion, when less frequent reports were common.
#' }
#' \item{`cases_tot`}{
#'   An integer vector. `cases_tot[i]` is the number of cases confirmed
#'   up to the end of `Date[i]` in `province[i]`.
#' }
#' \item{`cases_new`}{
#'   An integer vector. `cases_new[i]` is the number of cases confirmed
#'   from the end of `Date[i-1]` to the end of `Date[i]` in `province[i]`.
#' }
#' }
#'
#' @source
#' Raw data were downloaded from Michael Li's public
#' [Github repository](https://github.com/wzmli).
#' Up-to-date data can be downloaded
#' [here](https://wzmli.github.io/COVID19-Canada/COVID19_Canada.csv).
#'
#' @usage data(canadacovid)
#' @examples
#' data(canadacovid)
#' subset(canadacovid, province == "ON", -province)
#'
#' @name canadacovid
"canadacovid"

#' Husting wills
#'
#' Counts of wills probated in the Court of Husting
#' during four plague epidemics in 14th century London.
#'
#' @format
#' A data frame with 723 rows and 4 variables:
#' \describe{
#' \item{`Date`}{
#'   A Date vector. Spacing varies, as only nonzero counts of wills
#'   are included.
#' }
#' \item{`wills`}{
#'   An integer vector. `wills[i]` is the number of wills written
#'   on `Date[i]`.
#' }
#' \item{`outbreak`}{
#'   An ordered factor, splitting the time series by plague outbreak.
#'   Levels indicate the years in which outbreaks took place:
#'   1348, 1361, 1368, and 1375.
#' }
#' \item{`severity`}{
#'   An ordered factor, classifying outbreaks as `"minor"` or `"major"`.
#'   All 14th century outbreaks are classified as major, as the data
#'   are too sparse to distinguish between minor and major outbreaks,
#'   in contrast with wills probated in the
#'   [Prerogative Court of Canterbury][canterbury].
#' }
#' }
#'
#' @source
#' These data were transcribed from
#' \insertCite{Shar89;textual}{epigrowthfit}.
#'
#' @references
#' \insertRef{Shar89}{epigrowthfit}
#'
#' @usage data(husting)
#' @examples
#' data(husting)
#' subset(husting, outbreak == 1375, -outbreak)
#' @name husting
"husting"

#' Canterbury wills
#'
#' Counts of wills probated in the Prerogative Court of Canterbury
#' during 24 plague epidemics in 16th and 17th century London.
#'
#' @format
#' A data frame with 32681 rows and 4 variables:
#' \describe{
#' \item{`Date`}{
#'   A Date vector. Spacing varies, as only nonzero counts of wills
#'   are included.
#' }
#' \item{`wills`}{
#'   An integer vector. `wills[i]` is the number of wills written
#'   on `Date[i]`.
#' }
#' \item{`outbreak`}{
#'   An ordered factor, splitting the time series by plague outbreak.
#'   Levels indicate the years in which outbreaks took place:
#'   1563, 1578, ..., 1647, and 1665 (the Great Plague of London).
#' }
#' \item{`severity`}{
#'   An ordered factor, classifying outbreaks as `"minor"` or `"major"`.
#'   An outbreak is classified as major if and only if plague deaths
#'   per week per 1000 individuals exceeded 5 at least once.
#' }
#' }
#'
#' @source
#' These data were retrieved from the National Archives (UK) in 2018
#' using this
#' [research guide](https://www.nationalarchives.gov.uk/help-with-your-research/research-guides/wills-1384-1858/).
#'
#' @usage data(canterbury)
#' @examples
#' data(canterbury)
#' subset(canterbury, outbreak == 1665, -outbreak) # Great Plague of London
#' @name canterbury
"canterbury"

#' Burials in London parish registers
#'
#' Weekly counts of burials listed in extant parish registers
#' during 24 plague epidemics in 16th and 17th century London.
#'
#' @format
#' A data frame with 14198 rows and 4 variables:
#' \describe{
#' \item{`Date`}{
#'   A Date vector. Weekly, except for 340 instances of eight day
#'   spacing.
#' }
#' \item{`burials`}{
#'   An integer vector. `burials[i]` is the number of burials
#'   registered from the end of `Date[i-1]` to the end of `Date[i]`.
#' }
#' \item{`outbreak`}{
#'   An ordered factor, splitting the time series by plague outbreak.
#'   Levels indicate the years in which outbreaks took place:
#'   1563, 1578, ..., 1647, and 1665 (the Great Plague of London).
#' }
#' \item{`severity`}{
#'   An ordered factor, classifying outbreaks as `"minor"` or `"major"`.
#'   An outbreak is classified as major if and only if plague deaths
#'   per week per 1000 individuals exceeded 5 at least once.
#' }
#' }
#'
#' @source
#' These data were retrieved by
#' \insertCite{Cumm+16;textual}{epigrowthfit}
#' with permission from
#' [Ancestry.com](https://www.ancestry.com/).
#' They have been made publicly available
#' [here](http://www.neilcummins.com).
#'
#' @references
#' \insertRef{Cumm+16}{epigrowthfit}
#'
#' @usage data(londonparishes)
#' @examples
#' data(londonparishes)
#' subset(londonparishes, outbreak == 1665, -outbreak) # Great Plague of London
#' @name londonparishes
"londonparishes"

#' Plague deaths in the London Bills of Mortality
#'
#' Weekly counts of plague deaths recorded in the London Bills of
#' Mortality during 24 plague epidemics in 16th and 17th century London.
#'
#' @format
#' A data frame with 2041 rows and 6 variables:
#' \describe{
#' \item{`Date`}{
#'   A Date vector. Weekly within each level of `outbreak`,
#'   except for 78 instances of spacing not equal to 7 days.
#' }
#' \item{`deaths_all_causes`}{
#'   An integer vector. `deaths_all_causes[i]` is the number of deaths
#'   due to all causes from the end of `Date[i-1]` to the end of `Date[i]`.
#' }
#' \item{`deaths_plague`}{
#'   An integer vector. `deaths_plague[i]` is the number of deaths
#'   due to plague from the end of `Date[i-1]` to the end of `Date[i]`.
#' }
#' \item{`population`}{
#'   An integer vector. Estimated London population size.
#' }
#' \item{`outbreak`}{
#'   An ordered factor, splitting the time series by plague outbreak.
#'   Levels indicate the years in which outbreaks took place:
#'   1563, 1578, ..., 1647, and 1665 (the Great Plague of London).
#' }
#' \item{`severity`}{
#'   An ordered factor, classifying outbreaks as `"minor"` or `"major"`.
#'   An outbreak is classified as major if and only if plague deaths
#'   per week per 1000 individuals exceeded 5 at least once.
#' }
#' }
#'
#' @usage data(londonbills)
#' @examples
#' data(londonbills)
#' subset(londonbills, outbreak == 1665, -outbreak) # Great Plague of London
#'
#' @name londonbills
"londonbills"

#' Latent period distribution (pneumonic plague)
#'
#' Empirical distribution of the latent period of pneumonic plague.
#'
#' @format
#' A data frame with 12 rows and 2 variables:
#' \describe{
#' \item{`days`}{Latent period in days, from 1 to 12.}
#' \item{`relfreq`}{Observed relative frequency out of 224 cases.}
#' }
#'
#' @source
#' Obtained from Figure 1A in
#' \insertCite{GaniLeac04;textual}{epigrowthfit}.
#'
#' @references
#' \insertRef{GaniLeac04}{epigrowthfit}
#'
#' @usage data(plague_latent_period)
#' @examples
#' data(plague_latent_period)
#' plot(relfreq ~ days, data = plague_latent_period)
#'
#' @name plague_latent_period
"plague_latent_period"

#' Infectious period distribution (pneumonic plague)
#'
#' Empirical distribution of the infectious period of pneumonic plague.
#'
#' @format
#' A data frame with 12 rows and 2 variables:
#' \describe{
#' \item{`days`}{Infectious period in days, from 1 to 12.}
#' \item{`relfreq`}{Observed relative frequency out of 225 cases.}
#' }
#'
#' @source
#' Obtained from Figure 1B in
#' \insertCite{GaniLeac04;textual}{epigrowthfit}.
#'
#' @references
#' \insertRef{GaniLeac04}{epigrowthfit}
#'
#' @usage data(plague_infectious_period)
#' @examples
#' data(plague_infectious_period)
#' plot(relfreq ~ days, data = plague_infectious_period)
#'
#' @name plague_infectious_period
"plague_infectious_period"

#' Generation interval distribution (COVID-19)
#'
#' Gamma distribution of the COVID-19 generation interval fitted
#' to data from a cluster of 45 cases in Tianjin, China.
#'
#' @format
#' A list with 4 elements:
#' \describe{
#' \item{`shape`}{
#'   Estimated shape parameter. See [stats::dgamma()].
#' }
#' \item{`scale`}{
#'   Estimated scale parameter. See [stats::dgamma()].
#' }
#' \item{`breaks`}{
#'   An integer vector listing numbers of days. Equal to `0:20`.
#' }
#' \item{`probs`}{
#'   A numeric vector of length `length(breaks-1)`. `probs[i]`
#'   is the probability that the generation interval is between
#'   `breaks[i]` and `breaks[i+1]` days, conditional on `shape`
#'   and `scale`.
#'   Equal to `diff(pgamma(breaks, shape = shape, scale = scale))`.
#' }
#' }
#'
#' @source
#' `shape` and `scale` were computed from the mean and standard
#' deviation reported in \insertCite{Gany+20;textual}{epigrowthfit},
#' Table 4, Scenario 2.
#'
#' @references
#' \insertRef{Gany+20}{epigrowthfit}
#'
#' @usage data(plague_infectious_period)
#' @examples
#' data(covid_generation_interval)
#' x <- 10^seq(-2, log10(20), length.out = 150L)
#' fx <- with(covid_generation_interval,
#'   dgamma(x, shape = shape, scale = scale)
#' )
#' plot(x, fx, type = "l", ylim = c(0, max(fx[-1L])))
#'
#' @name covid_generation_interval
"covid_generation_interval"
