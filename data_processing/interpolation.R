## bilinear interpolation of ensemble forecasts to station locations
rm(list=ls())
data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"

## load station data (based on example.R)

library(RSQLite)
if (file.exists(paste0(data_dir, "observations/station_data"))) {
  db <- dbConnect(SQLite(), 'station_data')
}

st_temp <- dbGetQuery(db, 'SELECT * FROM temp') # get whole table 'temp'
str(st_temp)

st_meta <- dbGetQuery(db, 'SELECT * FROM meta') # get whole table 'meta'
str(st_meta)

dbDisconnect(db)


## load forecast data

library(ncdf4)

nc <- nc_open(paste0(data_dir, "forecasts/ecmwf_t2m_2007.nc"))
nc

fc_raw_temp <- ncvar_get(nc, "t2m")
str(fc_raw_temp) # [longitude,latitude,number,time]  

fc_raw_lat <- ncvar_get(nc, "latitude")
fc_raw_lon <- ncvar_get(nc, "longitude")
fc_raw_member <- ncvar_get(nc, "number")
fc_raw_time_unconverted <- ncvar_get(nc, "time")
fc_raw_validtime <- as.POSIXct(3600*fc_raw_time_unconverted, origin = "1900-01-01 00:00", tz = "UTC")
rm(fc_raw_time_unconverted)

# compute forecast initialization times (00 UTC forecasts, 36 and 48h ahead)
fc_raw_inittime <- fc_raw_validtime - 3600*c(36,48)

nc_close(nc)


## bilinear interpolation to station locations

## for 2007, generate required objects and perform bilinear interpolation

library(akima)
?bilinear

st_validtime <- as.POSIXct(paste0(st_temp$DATUM, sprintf("%02.f", st_temp$STUNDE)), 
                           tz = "UTC", format = "%Y%m%d%H")
st_stations <- st_temp$STATIONS_ID

# array with interpolated forecasts [stationID, member, validtime]
fc_interpolated_temp <- array(dim = c(length(st_meta$STATIONS_ID),
                                      length(fc_raw_member),
                                      length(fc_raw_validtime)))

# array with observed values [stationID, validtime]
obs_temp <- array(dim = c(length(st_meta$STATIONS_ID),
                          length(fc_raw_validtime)))

stationlist <- st_meta$STATIONS_ID

for(thisst in stationlist){
  thisst_pos <- which(stationlist == thisst)
  if(thisst_pos == 1){
    cat("------------------", "\n", 
        "   Year", 2007, "\n",
        "   .........", "\n")
  }
  if(thisst_pos %% 50 == 0){
    cat("station", thisst_pos, "of", length(stationlist), "starting at", paste(Sys.time()),"\n"); flush(stdout())
  }
  
  for(vtime in fc_raw_validtime){
    vtime_pos <- which(fc_raw_validtime == vtime)
    
    st_temp_pos <- which(st_validtime == vtime & st_stations == thisst)
    # st_temp contains no NA values, in case of missing observations, st_temp_pos is "integer(0)"
    if(length(st_temp_pos) == 0){
      obs_temp[thisst_pos,vtime_pos] <- NA
    } else{
      obs_temp[thisst_pos,vtime_pos] <- st_temp$TEMP[st_temp_pos]
    }
    
    for(mem in 1:50){
      thisst_lat <- st_meta$LATITUDE[which(st_meta$STATIONS_ID == thisst)]
      thisst_lon <- st_meta$LONGITUDE[which(st_meta$STATIONS_ID == thisst)]
      
      # find latitude and longitude of surrounding grib boxes
      lat_low <- max(fc_raw_lat[which(fc_raw_lat <= thisst_lat)])
      lat_high <- lat_low + 0.5
      lon_low <- max(fc_raw_lon[which(fc_raw_lon <= thisst_lon)])
      lon_high <- lon_low + 0.5
      x_coord <- c(lon_low, lon_high)
      y_coord <- c(lat_low, lat_high)
      
      # find corresponding positions in nc dimensions
      lat_pos <- c(which(fc_raw_lat == lat_low),which(fc_raw_lat == lat_high))
      lon_pos <- c(which(fc_raw_lon == lon_low),which(fc_raw_lon == lon_high))
      
      # bilinearly interpolate forecasts
      tmp <- bilinear(x = x_coord, y = y_coord, 
                      z = fc_raw_temp[lon_pos, lat_pos, mem, vtime_pos],
                      x0 = thisst_lon, y0 = thisst_lat)$z      
      # ... and convert to degrees Celsius  
      fc_interpolated_temp[thisst_pos, mem, vtime_pos] <- tmp - 273.15
    }
  }
}

## for years 2008-2016, load remaining nc files and append vectors / matrices / arrays
library(abind)

other_files <- paste0(data_dir, "forecasts/ecmwf_t2m_", 2008:2016, ".nc")

for(filename in other_files){
  cat("------------------", "\n", 
      "   Year", c(2008:2016)[which(other_files == filename)], "\n",
      "   .........", "\n")
  
  # load new nc file
  nc <- nc_open(filename)
  
  # t2m and "time" are changed compared to 2007
  fc_raw_temp <- ncvar_get(nc, "t2m")
  fc_raw_time_unconverted <- ncvar_get(nc, "time")
  fc_raw_validtime_append <- as.POSIXct(3600*fc_raw_time_unconverted, origin = "1900-01-01 00:00", tz = "UTC")
  rm(fc_raw_time_unconverted)
  
  nc_close(nc)
  
  # array with interpolated forecasts [stationID, member, validtime], to be appended later
  fc_interpolated_temp_append <- array(dim = c(length(st_meta$STATIONS_ID),
                                               length(fc_raw_member),
                                               length(fc_raw_validtime_append)))
  
  # array with observed values [stationID, validtime], to be appended later
  obs_temp_append <- array(dim = c(length(st_meta$STATIONS_ID),
                                   length(fc_raw_validtime_append)))
  
  
  for(thisst in stationlist){
    thisst_pos <- which(stationlist == thisst)
    if(thisst_pos %% 50 == 0){
      cat("station", thisst_pos, "of", length(stationlist), "starting at", paste(Sys.time()),"\n"); flush(stdout())
    }
    
    for(vtime in fc_raw_validtime_append){
      # last observation is recorded prior to last forecast valid time
      if(vtime > "2016-12-31"){
        next
      }
      vtime_pos <- which(fc_raw_validtime_append == vtime)
      
      st_temp_pos <- which(st_validtime == vtime & st_stations == thisst)
      # st_temp contains no NA values, in case of missing observations, st_temp_pos is "integer(0)"
      if(length(st_temp_pos) == 0){
        obs_temp_append[thisst_pos,vtime_pos] <- NA
      } else{
        obs_temp_append[thisst_pos,vtime_pos] <- st_temp$TEMP[st_temp_pos]
      }
      
      for(mem in 1:50){
        thisst_lat <- st_meta$LATITUDE[which(st_meta$STATIONS_ID == thisst)]
        thisst_lon <- st_meta$LONGITUDE[which(st_meta$STATIONS_ID == thisst)]
        
        # find latitude and longitude of surrounding grib boxes
        lat_low <- max(fc_raw_lat[which(fc_raw_lat <= thisst_lat)])
        lat_high <- lat_low + 0.5
        lon_low <- max(fc_raw_lon[which(fc_raw_lon <= thisst_lon)])
        lon_high <- lon_low + 0.5
        x_coord <- c(lon_low, lon_high)
        y_coord <- c(lat_low, lat_high)
        
        # find corresponding positions in nc dimensions
        lat_pos <- c(which(fc_raw_lat == lat_low),which(fc_raw_lat == lat_high))
        lon_pos <- c(which(fc_raw_lon == lon_low),which(fc_raw_lon == lon_high))
        
        # bilinearly interpolate forecasts
        tmp <- bilinear(x = x_coord, y = y_coord, 
                        z = fc_raw_temp[lon_pos, lat_pos, mem, vtime_pos],
                        x0 = thisst_lon, y0 = thisst_lat)$z      
        # ... and convert to degrees Celsius  
        fc_interpolated_temp_append[thisst_pos, mem, vtime_pos] <- tmp - 273.15
      }
    }
  }
  
  # append new objects to existing ones
  obs_temp <- abind(obs_temp, obs_temp_append, along = 2)
  fc_interpolated_temp <- abind(fc_interpolated_temp,
                                fc_interpolated_temp_append,
                                along = 3)
  fc_raw_validtime <- c(fc_raw_validtime, fc_raw_validtime_append)
  rm(obs_temp_append, fc_interpolated_temp_append, fc_raw_validtime_append)
}


## ------------------------------------------ ##

## write to netCDF file

# define dimensions
stationdim <- ncdim_def("station", "station_ID", as.integer(stationlist))
memberdim <- ncdim_def("member", "member_number", as.integer(1:50))
timedim <- ncdim_def(name = "time", vals = as.integer(fc_raw_validtime),
                     units = "seconds since 1970-01-01 00:00 UTC",
                     longname = "valid time of forecasts and observations, UTC")
# as.POSIXct(as.integer(fc_raw_validtime), tz = "UTC", origin = "1970-01-01 00:00")

# define variables
fillvalue <- NA

dlname <- "interpolated t2m ensemble forecast"
t2m_fc_def <- ncvar_def(name = "t2m_fc", units = "deg_C",
                        dim = list(stationdim,memberdim,timedim),
                        missval = fillvalue, longname = dlname,
                        prec="single")

dlname <- "t2m station observation"
t2m_obs_def <- ncvar_def(name = "t2m_obs", units = "deg_C",
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
ncfile_name <- paste0(data_dir,"data_interpolated.nc")
ncout <- nc_create(ncfile_name,
                   list(t2m_fc_def, t2m_obs_def, alt_def, lat_def, lon_def, id_def, location_def),
                   force_v4=T)

# put variables
ncvar_put(nc = ncout, varid = t2m_fc_def, vals = fc_interpolated_temp)
ncvar_put(nc = ncout, varid = t2m_obs_def, vals = obs_temp)
ncvar_put(nc = ncout, varid = alt_def, vals = st_meta$ALTITUDE)
ncvar_put(nc = ncout, varid = lat_def, vals = st_meta$LATITUDE)
ncvar_put(nc = ncout, varid = lon_def, vals = st_meta$LONGITUDE)
ncvar_put(nc = ncout, varid = id_def, vals = st_meta$STATIONS_ID)
ncvar_put(nc = ncout, varid = location_def, vals = st_meta$LOCATION)

nc_close(ncout)