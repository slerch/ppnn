library(RSQLite)

# setwd()

if (file.exists('station_data')) {
  db <- dbConnect(SQLite(), 'station_data')
}

dbListTables(db) # 2 tables: 'temp' and 'meta'
dbListFields(db, "temp") # 4 variables
dbListFields(db, "meta") # 5 variables

temp <-
  dbGetQuery(db, 'SELECT * FROM temp') # get whole table 'temp'
str(temp)

meta <-
  dbGetQuery(db, 'SELECT * FROM meta') # get whole table 'meta'
str(meta)

# get Potsdam data
stationID <- dbGetQuery(db,
                        'SELECT STATIONS_ID
                        FROM meta
                        WHERE LOCATION == "Potsdam"')
str(stationID)

if (identical(nrow(stationID), 1L)) {
  data <- dbGetQuery(db,
                     sprintf(
                       'SELECT *
                       FROM temp
                       WHERE STATIONS_ID == %i',
                       stationID$STATIONS_ID
                     ))
}
str(data)

dbDisconnect(db)
