## bilinear interpolation of ensemble forecasts to station locations
## code adapted for wind speed data

rm(list=ls())
data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"

## load station data (based on example.R)

library(RSQLite)
if (file.exists(paste0(data_dir, "observations/station_data_wind"))) {
  db <- dbConnect(SQLite(), paste0(data_dir, "observations/station_data_wind"))
}

st_wind <- dbGetQuery(db, 'SELECT * FROM wind') # get whole table 'wind'
str(st_wind)

st_meta <- dbGetQuery(db, 'SELECT * FROM meta') # get whole table 'meta'
str(st_meta)

dbDisconnect(db)

## load forecast data (contained in another file compared to temperature data, originally downloaded as auxiliary data!)

library(ncdf4)

nc <- nc_open(paste0(data_dir, "forecasts/auxiliary/processed/ecmwf_aux_surface_more_2007.nc"))
nc

fc_raw_u10 <- ncvar_get(nc, "u10")
str(fc_raw_u10) # [longitude,latitude,number,time]  
fc_raw_v10 <- ncvar_get(nc, "v10")
str(fc_raw_v10) # [longitude,latitude,number,time]  

# compute wind speed from u and v component
fc_raw_ws <- sqrt(fc_raw_u10^2 + fc_raw_v10^2)
str(fc_raw_ws)
rm(fc_raw_u10, fc_raw_v10)

fc_raw_lat <- ncvar_get(nc, "latitude")
fc_raw_lon <- ncvar_get(nc, "longitude")
fc_raw_member <- ncvar_get(nc, "number")
fc_raw_time_unconverted <- ncvar_get(nc, "time")
fc_raw_validtime <- as.POSIXct(3600*fc_raw_time_unconverted, origin = "1900-01-01 00:00", tz = "UTC")
rm(fc_raw_time_unconverted)

# compute forecast initialization times (00 UTC forecasts, 36 and 48h ahead)
fc_raw_inittime <- fc_raw_validtime - 3600*c(36,48)

nc_close(nc)

## interpolation to station locations
## changed compared to temperature forecasts:
## use data from closest grid point instead of bilinear interpolation
## (probably makes more sense for wind speed)

st_validtime <- as.POSIXct(paste0(st_wind$DATUM, sprintf("%02.f", st_wind$STUNDE)), 
                           tz = "UTC", format = "%Y%m%d%H")
st_stations <- st_wind$STATIONS_ID

# array with interpolated forecasts [stationID, member, validtime]
fc_interpolated_wind <- array(dim = c(length(st_meta$STATIONS_ID),
                                      length(fc_raw_member),
                                      length(fc_raw_validtime)))

# array with observed values [stationID, validtime]
obs_wind <- array(dim = c(length(st_meta$STATIONS_ID),
                          length(fc_raw_validtime)))


stationlist <- st_meta$STATIONS_ID

for(thisst in stationlist){
  thisst_pos <- which(stationlist == thisst)
  if(thisst_pos == 1){
    cat("------------------", "\n", 
        "   Year", 2007, "\n",
        "   .........", "\n")
  }
  if(thisst_pos %% 10 == 0){
    cat("station", thisst_pos, "of", length(stationlist), "starting at", paste(Sys.time()),"\n"); flush(stdout())
  }

  # station positions
  thisst_pos <- which(stationlist == thisst)
  thisst_lat <- st_meta$LATITUDE[which(st_meta$STATIONS_ID == thisst)]
  thisst_lon <- st_meta$LONGITUDE[which(st_meta$STATIONS_ID == thisst)]
  
  # find latitude and longitude of surrounding grib boxes
  lat_pos_mindiff <- order(abs(fc_raw_lat - thisst_lat))[1]
  lon_pos_mindiff <- order(abs(fc_raw_lon - thisst_lon))[1]
  
  for(vtime in fc_raw_validtime){
    vtime_pos <- which(fc_raw_validtime == vtime)
    
    st_wind_pos <- which(st_validtime == vtime & st_stations == thisst)
    # st_temp contains no NA values, in case of missing observations, st_temp_pos is "integer(0)"
    if(length(st_wind_pos) == 0){
      obs_wind[thisst_pos, vtime_pos] <- NA
    } else{
      obs_wind[thisst_pos, vtime_pos] <- st_wind$WIND[st_wind_pos]
    }
    
    # choose forecasts from nearest grid point
    fc_interpolated_wind[thisst_pos, , vtime_pos] <- fc_raw_ws[lon_pos_mindiff, lat_pos_mindiff, , vtime_pos]
  }
}


## for years 2008-2016, load remaining nc files and append vectors / matrices / arrays
library(abind)

other_files <- paste0(data_dir, "forecasts/auxiliary/processed/ecmwf_aux_surface_more_", 2008:2016, ".nc")

for(filename in other_files){
  cat("------------------", "\n", 
      "   Year", c(2008:2016)[which(other_files == filename)], "\n",
      "   .........", "\n")
  
  # load new nc file
  nc <- nc_open(filename)
  
  fc_raw_u10 <- ncvar_get(nc, "u10")
  # str(fc_raw_u10) # [longitude,latitude,number,time]  
  fc_raw_v10 <- ncvar_get(nc, "v10")
  # str(fc_raw_v10) # [longitude,latitude,number,time]  
  
  # compute wind speed from u and v component
  fc_raw_ws <- sqrt(fc_raw_u10^2 + fc_raw_v10^2)
  # str(fc_raw_ws)
  rm(fc_raw_u10, fc_raw_v10)
  
  fc_raw_time_unconverted <- ncvar_get(nc, "time")
  fc_raw_validtime_append <- as.POSIXct(3600*fc_raw_time_unconverted, origin = "1900-01-01 00:00", tz = "UTC")
  rm(fc_raw_time_unconverted)
  
  nc_close(nc)
  
  # array with interpolated forecasts [stationID, member, validtime], to be appended later
  fc_interpolated_wind_append <- array(dim = c(length(st_meta$STATIONS_ID),
                                               length(fc_raw_member),
                                               length(fc_raw_validtime_append)))
  
  # array with observed values [stationID, validtime], to be appended later
  obs_wind_append <- array(dim = c(length(st_meta$STATIONS_ID),
                                   length(fc_raw_validtime_append)))

  for(thisst in stationlist){
    
    thisst_pos <- which(stationlist == thisst)
    thisst_lat <- st_meta$LATITUDE[which(st_meta$STATIONS_ID == thisst)]
    thisst_lon <- st_meta$LONGITUDE[which(st_meta$STATIONS_ID == thisst)]
    
    # find latitude and longitude of surrounding grib boxes
    lat_pos_mindiff <- order(abs(fc_raw_lat - thisst_lat))[1]
    lon_pos_mindiff <- order(abs(fc_raw_lon - thisst_lon))[1]
    
    if(thisst_pos %% 50 == 0){
      cat("station", thisst_pos, "of", length(stationlist), "starting at", paste(Sys.time()),"\n"); flush(stdout())
    }
    
    for(vtime in fc_raw_validtime_append){
      # last observation is recorded prior to last forecast valid time
      # if(vtime > "2016-12-31"){
      #   next
      # }
      vtime_pos <- which(fc_raw_validtime_append == vtime)
      
      st_wind_pos <- which(st_validtime == vtime & st_stations == thisst)
      # st_temp contains no NA values, in case of missing observations, st_temp_pos is "integer(0)"
      if(length(st_wind_pos) == 0){
        obs_wind_append[thisst_pos,vtime_pos] <- NA
      } else{
        obs_wind_append[thisst_pos,vtime_pos] <- st_wind$WIND[st_wind_pos]
      }
      
      # choose forecasts from nearest grid point
      fc_interpolated_wind_append[thisst_pos, , vtime_pos] <- fc_raw_ws[lon_pos_mindiff, lat_pos_mindiff, , vtime_pos]
    }
  }
  
  # append new objects to existing ones
  fc_interpolated_wind <- abind(fc_interpolated_wind,
                                fc_interpolated_wind_append,
                                along = 3)
  
  obs_wind <- abind(obs_wind, obs_wind_append, along = 2)
  
  fc_raw_validtime <- c(fc_raw_validtime, fc_raw_validtime_append)
  
  rm(fc_interpolated_wind_append, obs_wind_append, fc_raw_validtime_append)
}

## ------------------------------------------ ##

## write to netCDF file

## from here: not yet run

# define dimensions
stationdim <- ncdim_def("station", "station_ID", as.integer(stationlist))
memberdim <- ncdim_def("member", "member_number", as.integer(1:50))
timedim <- ncdim_def(name = "time", vals = as.integer(fc_raw_validtime),
                     units = "seconds since 1970-01-01 00:00 UTC",
                     longname = "valid time of forecasts and observations, UTC")
# as.POSIXct(as.integer(fc_raw_validtime), tz = "UTC", origin = "1970-01-01 00:00")

# define variables
fillvalue <- NA

dlname <- "interpolated 10m wind speed ensemble forecast"
ws_fc_def <- ncvar_def(name = "ws_fc", units = "m s^-1",
                        dim = list(stationdim,memberdim,timedim),
                        missval = fillvalue, longname = dlname,
                        prec="single")

dlname <- "10m wind speed station observation"
ws_obs_def <- ncvar_def(name = "ws_obs", units = "m s^-1",
                         dim = list(stationdim,timedim),
                         missval = fillvalue, longname = dlname,
                         prec="single")

dlname <- "altitude of station"
alt_def <- ncvar_def(name = "station_alt", units = "m",
                     dim = list(stationdim),
                     missval = fillvalue, longname = dlname,
                     prec = "single")

dlname <- "latitude of station"
lat_def <- ncvar_def(name = "station_lat", units = "degrees north",
                     dim = list(stationdim),
                     missval = fillvalue, longname = dlname,
                     prec = "single")

dlname <- "longitude of station"
lon_def <- ncvar_def(name = "station_lon", units = "degrees east",
                     dim = list(stationdim),
                     missval = fillvalue, longname = dlname,
                     prec = "single")

dlname <- "station ID"
id_def <- ncvar_def(name = "station_id", units = "",
                    dim = list(stationdim),
                    missval = fillvalue, longname = dlname,
                    prec = "single")

# character variables (location names) require special attention
dimnchar <- ncdim_def("nchar", "", 1:36, create_dimvar=FALSE,
                      longname = "number of characters for locations")
location_def <- ncvar_def("station_loc", "", list(dimnchar,stationdim),
                          prec = "char", longname = "location of station")

## create nc file
ncfile_name <- paste0(data_dir, "/wind/", "data_interpolated_wind.nc")
ncout <- nc_create(ncfile_name,
                   list(ws_fc_def, ws_obs_def, alt_def, lat_def, lon_def, id_def, location_def),
                   force_v4=T)

# put variables
ncvar_put(nc = ncout, varid = ws_fc_def, vals = fc_interpolated_wind)
ncvar_put(nc = ncout, varid = ws_obs_def, vals = obs_wind)
ncvar_put(nc = ncout, varid = alt_def, vals = st_meta$ALTITUDE)
ncvar_put(nc = ncout, varid = lat_def, vals = st_meta$LATITUDE)
ncvar_put(nc = ncout, varid = lon_def, vals = st_meta$LONGITUDE)
ncvar_put(nc = ncout, varid = id_def, vals = st_meta$STATIONS_ID)
ncvar_put(nc = ncout, varid = location_def, vals = st_meta$LOCATION)

nc_close(ncout)