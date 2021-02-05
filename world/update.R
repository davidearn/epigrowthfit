# Retrieve full DELVE COVID-19 data set
loc <- paste0(
  "https://raw.githubusercontent.com/",
  "rs-delve/covid19_datasets/master/dataset/",
  "combined_dataset_latest.csv"
)
world_full <- read.csv(url(loc))

# Extract incidence
world <- world_full[c("DATE", "cases_new", "ISO", "country_name")]

# Clean up
names(world) <- c("date", "cases", "country_code", "country_name")
world$date <- as.Date(world$date)
world$cases <- as.integer(world$cases)
world$country_code <- factor(world$country_code)
world$country_name <- factor(world$country_name)
save(world, file = "world.RData")
