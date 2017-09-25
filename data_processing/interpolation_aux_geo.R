## interpolation of auxiliary model data to stations
##    by assigning values from closest grid points to stations

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

## load auxiliary model data

library(ncdf4)

nc <- nc_open(paste0(data_dir, "forecasts/auxiliary/ecmwf_aux_geo_2007.nc"))
nc

# contains lsm [units = (0 - 1)], orog [units = m]
# note: 3 dimensions (longitude, latitude, time)
# should not change over time -> time dimension not required

fc_raw_lsm <- ncvar_get(nc, "lsm")
str(fc_raw_lsm) # [longitude,latitude,number,time]  
fc_raw_orog <- ncvar_get(nc, "orog")

fc_raw_lsm_use <- fc_raw_lsm[,,1]
fc_raw_orog_use <- fc_raw_orog[,,1]

fc_raw_lat <- ncvar_get(nc, "latitude")
fc_raw_lon <- ncvar_get(nc, "longitude")

nc_close(nc)

## interpolation

# st_validtime <- as.POSIXct(paste0(st_temp$DATUM, sprintf("%02.f", st_temp$STUNDE)), 
#                            tz = "UTC", format = "%Y%m%d%H")
st_stations <- st_temp$STATIONS_ID

# array with interpolated forecasts [stationID, member, validtime]
fc_interpolated_lsm <- array(dim = c(length(st_meta$STATIONS_ID)))
fc_interpolated_orog <- fc_interpolated_lsm

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
  
  fc_interpolated_lsm[thisst_pos] <- fc_raw_lsm_use[lon_pos_mindiff, lat_pos_mindiff]
  fc_interpolated_orog[thisst_pos] <- fc_raw_orog_use[lon_pos_mindiff, lat_pos_mindiff]
}


## no need to go through years 2008-2016, lsm and orog don't change


## ------------------------------------------ ##

## write to netCDF file

# define dimensions
stationdim <- ncdim_def("station", "station_ID", as.integer(stationlist))

# define variables
fillvalue <- NA

dlname <- "interpolated land sea mask"
lsm_fc_def <- ncvar_def(name = "lsm", units = "(0 - 1)",
                        dim = list(stationdim),
                        missval = fillvalue, longname = dlname,
                        prec="single")

dlname <- "interpolated orography information"
orog_fc_def <- ncvar_def(name = "orog", units = "m",
                         dim = list(stationdim),
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
ncfile_name <- paste0(data_dir,"data_aux_geo_interpolated.nc")
ncout <- nc_create(ncfile_name,
                   list(lsm_fc_def, orog_fc_def, alt_def, lat_def, lon_def, id_def, location_def),
                   force_v4=T)

# put variables
ncvar_put(nc = ncout, varid = lsm_fc_def, vals = fc_interpolated_lsm)
ncvar_put(nc = ncout, varid = orog_fc_def, vals = fc_interpolated_orog)
ncvar_put(nc = ncout, varid = alt_def, vals = st_meta$ALTITUDE)
ncvar_put(nc = ncout, varid = lat_def, vals = st_meta$LATITUDE)
ncvar_put(nc = ncout, varid = lon_def, vals = st_meta$LONGITUDE)
ncvar_put(nc = ncout, varid = id_def, vals = st_meta$STATIONS_ID)
ncvar_put(nc = ncout, varid = location_def, vals = st_meta$LOCATION)

nc_close(ncout)

