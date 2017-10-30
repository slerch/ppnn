## interpolation of auxiliary model data to stations
##    by assigning values from closest grid points to stations

## split into two parts, otherwise R crashes due to lack of available memory
## this is part 2

rm(list=ls())
data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"

library(RSQLite)
if (file.exists(paste0(data_dir, "observations/station_data"))) {
  db <- dbConnect(SQLite(), paste0(data_dir, "observations/station_data"))
}

st_temp <- dbGetQuery(db, 'SELECT * FROM temp') # get whole table 'temp'
str(st_temp)

st_meta <- dbGetQuery(db, 'SELECT * FROM meta') # get whole table 'meta'
str(st_meta)

dbDisconnect(db)

## load (new set of) auxiliary model data

library(ncdf4)

nc <- nc_open(paste0(data_dir, "forecasts/auxiliary/ecmwf_aux_surface_more_2007.nc"))
nc

# variables:
#   d2m (2m Dew Point Temperature); units: K
#   ssr (Net short wave radiation flux - accumulated _ surface); units: W m-2
#   str (Net long wave radiation flux - accumulated _  surface); units: W m-2
#   sm (Soil Moisture); units: kg m**-3

fc_raw_d2m <- ncvar_get(nc, "d2m")
fc_raw_ssr <- ncvar_get(nc, "ssr")
fc_raw_str <- ncvar_get(nc, "str")
fc_raw_sm <- ncvar_get(nc, "sm")

fc_raw_lat <- ncvar_get(nc, "latitude")
fc_raw_lon <- ncvar_get(nc, "longitude")
fc_raw_member <- ncvar_get(nc, "number")
fc_raw_time_unconverted <- ncvar_get(nc, "time")
fc_raw_validtime <- as.POSIXct(3600*fc_raw_time_unconverted, origin = "1900-01-01 00:00", tz = "UTC")
rm(fc_raw_time_unconverted)

# compute forecast initialization times (00 UTC forecasts, 36 and 48h ahead)
fc_raw_inittime <- fc_raw_validtime - 3600*c(36,48)

nc_close(nc)

## interpolation

st_validtime <- as.POSIXct(paste0(st_temp$DATUM, sprintf("%02.f", st_temp$STUNDE)), 
                           tz = "UTC", format = "%Y%m%d%H")
st_stations <- st_temp$STATIONS_ID

# array with interpolated forecasts [stationID, member, validtime]
fc_interpolated_d2m <- array(dim = c(length(st_meta$STATIONS_ID),
                                      length(fc_raw_member),
                                      length(fc_raw_validtime)))
fc_interpolated_ssr <- fc_interpolated_str <- fc_interpolated_sm <- fc_interpolated_d2m


stationlist <- st_meta$STATIONS_ID


for(thisst in stationlist){
  
  thisst_pos <- which(stationlist == thisst)
  thisst_lat <- st_meta$LATITUDE[which(st_meta$STATIONS_ID == thisst)]
  thisst_lon <- st_meta$LONGITUDE[which(st_meta$STATIONS_ID == thisst)]
  
  # find latitude and longitude of surrounding grib boxes
  lat_pos_mindiff <- order(abs(fc_raw_lat - thisst_lat))[1]
  lon_pos_mindiff <- order(abs(fc_raw_lon - thisst_lon))[1]
  
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
    
    # choose forecasts from nearest grid point
    # d2m
    fc_interpolated_d2m[thisst_pos, , vtime_pos] <- fc_raw_d2m[lon_pos_mindiff, lat_pos_mindiff, , vtime_pos]
    # ssr
    fc_interpolated_ssr[thisst_pos, , vtime_pos] <- fc_raw_ssr[lon_pos_mindiff, lat_pos_mindiff, , vtime_pos]
    # str
    fc_interpolated_str[thisst_pos, , vtime_pos] <- fc_raw_str[lon_pos_mindiff, lat_pos_mindiff, , vtime_pos]
    # sm
    fc_interpolated_sm[thisst_pos, , vtime_pos] <- fc_raw_sm[lon_pos_mindiff, lat_pos_mindiff, , vtime_pos]
  }
}

## for years 2008-2016, load remaining nc files and append vectors / matrices / arrays
library(abind)

other_files <- paste0(data_dir, "forecasts/auxiliary/ecmwf_aux_surface_more_", 2008:2016, ".nc")

for(filename in other_files){
  cat("------------------", "\n", 
      "   Year", c(2008:2016)[which(other_files == filename)], "\n",
      "   .........", "\n")
  
  # load new nc file
  nc <- nc_open(filename)
  
  fc_raw_d2m <- ncvar_get(nc, "d2m")
  fc_raw_ssr <- ncvar_get(nc, "ssr")
  fc_raw_str <- ncvar_get(nc, "str")
  fc_raw_sm <- ncvar_get(nc, "sm")
  
  fc_raw_time_unconverted <- ncvar_get(nc, "time")
  fc_raw_validtime_append <- as.POSIXct(3600*fc_raw_time_unconverted, origin = "1900-01-01 00:00", tz = "UTC")
  rm(fc_raw_time_unconverted)
  
  nc_close(nc)
  
  # array with interpolated forecasts [stationID, member, validtime], to be appended later
  fc_interpolated_d2m_append <- array(dim = c(length(st_meta$STATIONS_ID),
                                               length(fc_raw_member),
                                               length(fc_raw_validtime_append)))
  fc_interpolated_ssr_append <- fc_interpolated_str_append <- fc_interpolated_sm_append <- fc_interpolated_d2m_append
  
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
      
      # choose forecasts from nearest grid point
      # d2m
      fc_interpolated_d2m_append[thisst_pos, , vtime_pos] <- fc_raw_d2m[lon_pos_mindiff, lat_pos_mindiff, , vtime_pos]
      # ssr
      fc_interpolated_ssr_append[thisst_pos, , vtime_pos] <- fc_raw_ssr[lon_pos_mindiff, lat_pos_mindiff, , vtime_pos]
      # str
      fc_interpolated_str_append[thisst_pos, , vtime_pos] <- fc_raw_str[lon_pos_mindiff, lat_pos_mindiff, , vtime_pos]
      # sm
      fc_interpolated_sm_append[thisst_pos, , vtime_pos] <- fc_raw_sm[lon_pos_mindiff, lat_pos_mindiff, , vtime_pos]
      
    }
  }
  
  # append new objects to existing ones
  fc_interpolated_d2m <- abind(fc_interpolated_d2m,
                               fc_interpolated_d2m_append,
                               along = 3)
  fc_interpolated_ssr <- abind(fc_interpolated_ssr,
                               fc_interpolated_ssr_append,
                               along = 3)
  fc_interpolated_str <- abind(fc_interpolated_str,
                               fc_interpolated_str_append,
                               along = 3)
  fc_interpolated_sm <- abind(fc_interpolated_sm,
                              fc_interpolated_sm_append,
                              along = 3)
  
  fc_raw_validtime <- c(fc_raw_validtime, fc_raw_validtime_append)
  
  rm(fc_interpolated_d2m_append, fc_interpolated_ssr_append, fc_interpolated_str_append, fc_interpolated_sm_append)
}

## ------------------------------------------ ##

## write to netCDF file

# define dimensions
stationdim <- ncdim_def("station", "station_ID", as.integer(stationlist))
memberdim <- ncdim_def("member", "member_number", as.integer(1:50))
timedim <- ncdim_def(name = "time", vals = as.integer(fc_raw_validtime),
                     units = "seconds since 1970-01-01 00:00 UTC",
                     longname = "valid time of forecasts and observations, UTC")

# define variables
fillvalue <- NA

dlname <- "interpolated 2m dew point temperature ensemble forecast"
d2m_fc_def <- ncvar_def(name = "d2m_fc", units = "K",
                         dim = list(stationdim,memberdim,timedim),
                         missval = fillvalue, longname = dlname,
                         prec="single")

dlname <- "interpolated net short wave radiation flux ensemble forecast"
ssr_fc_def <- ncvar_def(name = "ssr_fc", units = "W m^-2",
                         dim = list(stationdim,memberdim,timedim),
                         missval = fillvalue, longname = dlname,
                         prec="single")

dlname <- "interpolated net long wave radiation flux ensemble forecast"
str_fc_def <- ncvar_def(name = "str_fc", units = "W m^-2",
                        dim = list(stationdim,memberdim,timedim),
                        missval = fillvalue, longname = dlname,
                        prec="single")

dlname <- "interpolated soil moisture ensemble forecast"
sm_fc_def <- ncvar_def(name = "sm_fc", units = "kg m^-3",
                        dim = list(stationdim,memberdim,timedim),
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
ncfile_name <- paste0(data_dir,"data_aux_surface_more_interpolated_part2.nc")
ncout <- nc_create(ncfile_name,
                   list(ssr_fc_def, str_fc_def, d2m_fc_def, sm_fc_def, 
                        alt_def, lat_def, lon_def, id_def, location_def),
                   force_v4=T)

# put variables
ncvar_put(nc = ncout, varid = ssr_fc_def, vals = fc_interpolated_ssr)
ncvar_put(nc = ncout, varid = str_fc_def, vals = fc_interpolated_str)
ncvar_put(nc = ncout, varid = d2m_fc_def, vals = fc_interpolated_d2m)
ncvar_put(nc = ncout, varid = sm_fc_def, vals = fc_interpolated_sm)

ncvar_put(nc = ncout, varid = alt_def, vals = st_meta$ALTITUDE)
ncvar_put(nc = ncout, varid = lat_def, vals = st_meta$LATITUDE)
ncvar_put(nc = ncout, varid = lon_def, vals = st_meta$LONGITUDE)
ncvar_put(nc = ncout, varid = id_def, vals = st_meta$STATIONS_ID)
ncvar_put(nc = ncout, varid = location_def, vals = st_meta$LOCATION)

nc_close(ncout)
