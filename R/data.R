#' Plague data
#'
#' @name plague
#' @aliases data deaths mortality
#' @format A data frame with 1429 observations and 12 variables:
#' \describe{
#' \item{year}{year}
#' \item{month}{month}
#' \item{day}{day}
#' \item{time}{decimal date (year + fraction-of-year)}
#' \item{plague.deaths}{recorded plague deaths}
#' \item{all.cause.deaths}{recorded deaths from all causes}
#' \item{aggregation}{time unit (e.g., monthly, weekly)}
#' \item{place}{where the epidemic occurred}
#' \item{outbreak.year}{the year in which the epidemic took off (not always the year of the first report)}
#' \item{type}{what was recorded (e.g., testaments, deaths)}
#' \item{population}{estimate of population size, if available}
#' \item{severity}{minor vs major epidemic}
#' }
#'
#' @description
#' Population level data for plague epidemics in London, England, and Harbin, China.
#' The early data from the 14th century in London are monthly counts of extant
#' last wills and testaments, whereas all the later data come from published
#' weekly mortality records.
#'
#' These data are also available online from the International Infectious Disease Data Archive (IIDDA) at \href{http://iidda.mcmaster.ca}{http://iidda.mcmaster.ca}; see the IIDDA website for a list of all papers that have used these data.
#' @author
#' These data were compiled from Cohn (2003), Creighton (1965) and Dietz (2009).
#' @source
#' \itemize{
#' \item London, 1348-1349:  Cohn (2003, Figure 7.33, page 184).
#' \item London, 1361: \emph{ibid.}, Figure 7.34.
#' \item London, 1375: \emph{ibid.}, Figure 7.35.
#' \item London, 1563--1583: Creighton (1965)
#' \item London, 1593--1666: London Bills of Mortality
#' \item Harbin, 1911--1912: Dietz (2009)
#' }
#' @references
#' \insertRef{Cohn03}{epigrowthfitPNAS}
#'
#' \insertRef{Crei65}{epigrowthfitPNAS}
#'
#' \insertRef{Diet09}{epigrowthfitPNAS}
#' @examples
#' ## Look at the first few records in the data frame
#' head(plague)
"plague"

#' Estimated infectious period for pneumonic plague
#'
#' @name infectious.period
#' @format A data frame containing columns:
#' \describe{
#' \item{\code{days}}{days since beginning of infectiousness}
#' \item{\code{frequency}}{frequency distribution of the length of the infectious period}
#' \item{\code{prediction}}{fitted lognormal distribution}
#' }
#' @source
#' These data come from Figure 1B of Gani and Leach (2004).
#' @references
#' \insertRef{GaniLeac04}{epigrowthfitPNAS}
"infectious.period"

#' Estimated latent period for pneumonic plague
#'
#' @name latent.period
#' @format A data frame containing columns:
#' \describe{
#' \item{\code{days}}{days since infection}
#' \item{\code{frequency}}{frequency distribution of the length of the latent period}
#' \item{\code{prediction}}{fitted lognormal distribution}
#' }
#' @source
#' These data come from Figure 1A of Gani and Leach (2004).
#' @references
#' \insertRef{GaniLeac04}{epigrowthfitPNAS}
"latent.period"

#' Epidemic year definitions
#'
#' @name epidemic_defs
#' @format A data frame containing columns:
#' \describe{
#' \item{\code{outbreak.year}}{numeric outbreak year}
#' \item{\code{severity}}{"major" or "minor"}
#' \item{\code{start}}{starting time (julian date)}
#' \item{\code{end}}{ending time (julian date)}
#' }
"epidemic_defs"

#' Canada COVID-19 cases, 2020
#'
#' @name canada_covid
#' @format A tibble containing columns:
#' \describe{
#' \item{\code{Date}}{reporting date}
#' \item{\code{Province}}{province}
#' \time{\code{newConfirmations}}{number of new confirmed cases}
#' }
"canada_covid"
