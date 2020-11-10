url <- "https://wzmli.github.io/COVID19-Canada/COVID19_Canada.csv"
filename <- "../canadacovid_raw.csv"
download.file(url, filename)
canadacovid_raw <- read.csv(filename, header = TRUE)
save(canadacovid_raw, file = sub(".csv", ".RData", filename))
file.remove(filename)
