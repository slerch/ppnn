## based on code by Alexander Jordan, adapted for wind speed data

## note the terms of use: ftp://ftp-cdc.dwd.de/pub/CDC/Terms_of_use.pdf

library(RCurl)
library(magrittr)
library(stringr)

#### paths ####
myPath <- "/path/to/data/"

# download folder
dlPath <-
  paste0(myPath, "data_files_wind/download/")
# data folder
dataPath <- 
  paste0(myPath, "data_files_wind/")
# download path
dwdPath <-
  "ftp://ftp-cdc.dwd.de/pub/CDC/observations_germany/climate/hourly/wind/historical/"

# create folders if necessary
if (!dir.exists(dlPath)) dir.create(dlPath, recursive = TRUE)


#### download files ####

fileNames <- dwdPath %>%
  getURL(dirlistonly = TRUE) %>%
  strsplit("\n") %>%
  .[[1]]

lapply(fileNames, function(name) {
  url <- paste0(dwdPath, name)
  desturl <- paste0(dlPath, name)
  if (!file.exists(desturl)) download.file(url, desturl, mode = "wb")
})

#### unzip/copy files ####

dirList <- dir(dlPath)
zipped <- grepl(".zip", dirList)
zipList <- str_c(dlPath, dirList[zipped])
lapply(zipList, unzip, exdir = dataPath)
copyList <- dirList[!zipped]
file.copy(str_c(dlPath, copyList),
          str_c(dataPath, "/", copyList),
          overwrite = TRUE)