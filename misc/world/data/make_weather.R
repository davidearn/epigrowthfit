(dir_weather <- Sys.getenv("DIR_WEATHER"))
(path_utils  <- Sys.getenv("PATH_UTILS"))
(path_coords <- Sys.getenv("PATH_COORDS"))
(path_rds    <- Sys.getenv("PATH_RDS"))

options(timeout = 600, con_total = 100L, con_host = 6L, warn = 1L, mc.cores = 4L)
set.seed(960850L, kind = "L'Ecuyer-CMRG")

utils <- new.env()
source(path_utils, local = utils)

coords <- readRDS(path_coords)
coords <- coords[c("country_iso_alpha3", "longitude", "latitude", "population")]
names(coords)[1L] <- "country"

weather_long <- evalq(envir = utils, {
  update_weather(
    path = dir_weather,
    data = coords,
    n = 1000L,
    scale = 3
  )
  shape_weather(dir_weather)
})
names(weather_long)[match("country", names(weather_long), 0L)] <- "country_iso_alpha3"
levels(weather_long[["varid"]]) <- paste0("weather_", levels(weather_long[["varid"]]))

weather_wide <- reshape(
  data = weather_long,
  direction = "wide",
  varying = levels(weather_long[["varid"]]),
  v.names = "value",
  timevar = "varid",
  idvar = c("country_iso_alpha3", "date")
)
attr(weather_wide, "reshapeWide") <- NULL

saveRDS(weather_wide, file = path_rds)
str(weather_wide)
