## based on code by Alexander Jordan, adapted for wind speed data

rm(list=ls())

library(RSQLite)
library(stringr)
library(readr)
library(dplyr)
library(lubridate)

#### paths ####
myPath <- "/path/to/data/"

# data folder
dataPath <- 
  paste0(myPath, "data_files_wind/")
# database file
dbPath <-
  paste0(myPath, "data_files/station_data_wind")

#### build database ####

# wind filename structure "produkt_ff_stunde_'START'_'END'_'ID'.txt"
# START 19:26
# END 28:35
# ID 37:41
getEND <- function(x) str_sub(x, 28, 35)
getID <- function(x) str_sub(x, 37, 41)

tempCsvList <- dir(dataPath, pattern = "produkt_ff_stunde_") %>%
  (function(x) x[as.integer(getEND(x)) >= 20070101L]) %>%
  (function(x) x[order(as.integer(getID(x)))])

# build "wind" table
db <- dbConnect(SQLite(), dbPath)
lapply(seq_along(tempCsvList), function(i) {
  try({
    read_delim(
      str_c(dataPath, tempCsvList[i]),
      delim = ";",
      col_types = cols_only(
        STATIONS_ID = col_integer(),
        MESS_DATUM = col_datetime(format = "%Y%m%d%H"),
        F = col_number()
      ),
      trim_ws = TRUE
    ) %>%
      filter(year(MESS_DATUM) %in% 2007:2016) %>%
      filter(hour(MESS_DATUM) %in% c(0, 12)) %>%
      filter(F != -999) %>%
      mutate(DATUM = as.integer(format(MESS_DATUM, format = "%Y%m%d"))) %>%
      mutate(STUNDE = as.integer(format(MESS_DATUM, format = "%H"))) %>%
      mutate(WIND = F) %>%
      mutate(MESS_DATUM = NULL) %>%
      mutate(F = NULL) %>%
      dbWriteTable(db,
                   name = "wind",
                   value = .,
                   append = TRUE)
    
    message(sprintf("%d of %d", i, length(tempCsvList)))
  })
})
dbDisconnect(db)

# build "meta" table
db <- dbConnect(SQLite(), dbPath)
try({
  IDs <- dbGetQuery(db, 'SELECT STATIONS_ID FROM wind') %>%
    .$STATIONS_ID %>%
    unique
  
  read_fwf(
    str_c(dataPath, "FF_Stundenwerte_Beschreibung_Stationen.txt"),
    fwf_positions(
      c(1, 25, 40, 52, 62),
      c(5, 38, 50, 60, 101),
      c("STATIONS_ID", "ALTITUDE", "LATITUDE", "LONGITUDE", "LOCATION")
    ),
    col_types = cols_only(
      STATIONS_ID = "i",
      ALTITUDE = "n",
      LATITUDE = "n",
      LONGITUDE = "n",
      LOCATION = "c"
    ),
    skip = 2
  ) %>%
    filter(STATIONS_ID %in% IDs) %>%
    dbWriteTable(db,
                 name = "meta",
                 value = .,
                 append = TRUE)
})
dbDisconnect(db)