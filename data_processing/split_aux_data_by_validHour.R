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

library(ncdf4)

# ------------------------------------------ #

### pl 500 hPa data ###

nc <- nc_open(paste0(data_dir,"data_aux_pl500_interpolated.nc"))
nc

# u_pl500_fc, v_pl500_fc, gh_pl500_fc

u_pl500_fc <- ncvar_get(nc, "u_pl500_fc")
v_pl500_fc <- ncvar_get(nc, "v_pl500_fc")
gh_pl500_fc <- ncvar_get(nc, "gh_pl500_fc")

station_alt <- ncvar_get(nc, "station_alt")
station_lat <- ncvar_get(nc, "station_lat")
station_lon <- ncvar_get(nc, "station_lon")
station_id <- ncvar_get(nc, "station_id")
station_loc <- ncvar_get(nc, "station_loc")

time <- ncvar_get(nc, "time")
station <- ncvar_get(nc, "station")

nc_close(nc)

validtime <- as.POSIXct(time, tz = "UTC", origin = "1970-01-01 00:00")
validtime_hour <- format(validtime, "%H")

ind00 <- which(validtime_hour == "00")
ind12 <- which(validtime_hour == "12")

u_pl500_fc_00 <- u_pl500_fc[,,ind00]
u_pl500_fc_12 <- u_pl500_fc[,,ind12]

v_pl500_fc_00 <- v_pl500_fc[,,ind00]
v_pl500_fc_12 <- v_pl500_fc[,,ind12]

gh_pl500_fc_00 <- gh_pl500_fc[,,ind00]
gh_pl500_fc_12 <- gh_pl500_fc[,,ind12]

time_00 <- time[ind00]
time_12 <- time[ind12]

rm(time, u_pl500_fc, v_pl500_fc, gh_pl500_fc)

## write to nc files

## 00 UTC file

# define dimensions
stationdim <- ncdim_def("station", "station_ID", station_id)
memberdim <- ncdim_def("member", "member_number", as.integer(1:50))
timedim <- ncdim_def(name = "time", vals = as.integer(time_00),
                     units = "seconds since 1970-01-01 00:00 UTC",
                     longname = "valid time of forecasts and observations, UTC")


# define variables
fillvalue <- NA

dlname <- "interpolated u wind component at pressure level 500 hPa ensemble forecast"
u_fc_def <- ncvar_def(name = "u_pl500_fc", units = "m s^-1",
                      dim = list(stationdim,memberdim,timedim),
                      missval = fillvalue, longname = dlname,
                      prec="single")

dlname <- "interpolated v wind component at pressure level 500 hPa ensemble forecast"
v_fc_def <- ncvar_def(name = "v_pl500_fc", units = "m s^-1",
                      dim = list(stationdim,memberdim,timedim),
                      missval = fillvalue, longname = dlname,
                      prec="single")

dlname <- "interpolated geopotential height at pressure level 500 hPa ensemble forecast"
gh_fc_def <- ncvar_def(name = "gh_pl500_fc", units = "gpm",
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
ncfile_name <- paste0(data_dir,"data_aux_pl500_interpolated_00UTC.nc")
ncout <- nc_create(ncfile_name,
                   list(u_fc_def, v_fc_def, gh_fc_def, alt_def, lat_def, lon_def, id_def, location_def),
                   force_v4=T)

# put variables
ncvar_put(nc = ncout, varid = u_fc_def, vals = u_pl500_fc_00)
ncvar_put(nc = ncout, varid = v_fc_def, vals = v_pl500_fc_00)
ncvar_put(nc = ncout, varid = gh_fc_def, vals = gh_pl500_fc_00)
ncvar_put(nc = ncout, varid = alt_def, vals = st_meta$ALTITUDE)
ncvar_put(nc = ncout, varid = lat_def, vals = st_meta$LATITUDE)
ncvar_put(nc = ncout, varid = lon_def, vals = st_meta$LONGITUDE)
ncvar_put(nc = ncout, varid = id_def, vals = st_meta$STATIONS_ID)
ncvar_put(nc = ncout, varid = location_def, vals = st_meta$LOCATION)

nc_close(ncout)


## 12 UTC file

# define dimensions
rm(timedim)
timedim <- ncdim_def(name = "time", vals = as.integer(time_12),
                     units = "seconds since 1970-01-01 00:00 UTC",
                     longname = "valid time of forecasts and observations, UTC")


# define variables
dlname <- "interpolated u wind component at pressure level 500 hPa ensemble forecast"
u_fc_def <- ncvar_def(name = "u_pl500_fc", units = "m s^-1",
                      dim = list(stationdim,memberdim,timedim),
                      missval = fillvalue, longname = dlname,
                      prec="single")

dlname <- "interpolated v wind component at pressure level 500 hPa ensemble forecast"
v_fc_def <- ncvar_def(name = "v_pl500_fc", units = "m s^-1",
                      dim = list(stationdim,memberdim,timedim),
                      missval = fillvalue, longname = dlname,
                      prec="single")

dlname <- "interpolated geopotential height at pressure level 500 hPa ensemble forecast"
gh_fc_def <- ncvar_def(name = "gh_pl500_fc", units = "gpm",
                       dim = list(stationdim,memberdim,timedim),
                       missval = fillvalue, longname = dlname,
                       prec="single")


## create nc file
ncfile_name <- paste0(data_dir,"data_aux_pl500_interpolated_12UTC.nc")
ncout <- nc_create(ncfile_name,
                   list(u_fc_def, v_fc_def, gh_fc_def, alt_def, lat_def, lon_def, id_def, location_def),
                   force_v4=T)

# put variables
ncvar_put(nc = ncout, varid = u_fc_def, vals = u_pl500_fc_12)
ncvar_put(nc = ncout, varid = v_fc_def, vals = v_pl500_fc_12)
ncvar_put(nc = ncout, varid = gh_fc_def, vals = gh_pl500_fc_12)
ncvar_put(nc = ncout, varid = alt_def, vals = st_meta$ALTITUDE)
ncvar_put(nc = ncout, varid = lat_def, vals = st_meta$LATITUDE)
ncvar_put(nc = ncout, varid = lon_def, vals = st_meta$LONGITUDE)
ncvar_put(nc = ncout, varid = id_def, vals = st_meta$STATIONS_ID)
ncvar_put(nc = ncout, varid = location_def, vals = st_meta$LOCATION)

nc_close(ncout)


# ------------------------------------------ #

### pl 850 hPa data ###

data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
nc <- nc_open(paste0(data_dir,"data_aux_pl850_interpolated.nc"))
nc

# u_pl850_fc, v_pl850_fc, q_pl850_fc

u_pl850_fc <- ncvar_get(nc, "u_pl850_fc")
v_pl850_fc <- ncvar_get(nc, "v_pl850_fc")
q_pl850_fc <- ncvar_get(nc, "q_pl850_fc")

station_alt <- ncvar_get(nc, "station_alt")
station_lat <- ncvar_get(nc, "station_lat")
station_lon <- ncvar_get(nc, "station_lon")
station_id <- ncvar_get(nc, "station_id")
station_loc <- ncvar_get(nc, "station_loc")

time <- ncvar_get(nc, "time")
station <- ncvar_get(nc, "station")

nc_close(nc)

validtime <- as.POSIXct(time, tz = "UTC", origin = "1970-01-01 00:00")
validtime_hour <- format(validtime, "%H")

ind00 <- which(validtime_hour == "00")
ind12 <- which(validtime_hour == "12")

u_pl850_fc_00 <- u_pl850_fc[,,ind00]
u_pl850_fc_12 <- u_pl850_fc[,,ind12]

v_pl850_fc_00 <- v_pl850_fc[,,ind00]
v_pl850_fc_12 <- v_pl850_fc[,,ind12]

q_pl850_fc_00 <- q_pl850_fc[,,ind00]
q_pl850_fc_12 <- q_pl850_fc[,,ind12]

time_00 <- time[ind00]
time_12 <- time[ind12]

rm(time, u_pl850_fc, v_pl850_fc, q_pl850_fc)

## write to nc files

## 00 UTC file

# define dimensions
stationdim <- ncdim_def("station", "station_ID", station_id)
memberdim <- ncdim_def("member", "member_number", as.integer(1:50))
timedim <- ncdim_def(name = "time", vals = as.integer(time_00),
                     units = "seconds since 1970-01-01 00:00 UTC",
                     longname = "valid time of forecasts and observations, UTC")

# define variables
fillvalue <- NA

dlname <- "interpolated u wind component at pressure level 850 hPa ensemble forecast"
u_fc_def <- ncvar_def(name = "u_pl850_fc", units = "m s^-1",
                      dim = list(stationdim,memberdim,timedim),
                      missval = fillvalue, longname = dlname,
                      prec="single")

dlname <- "interpolated v wind component at pressure level 850 hPa ensemble forecast"
v_fc_def <- ncvar_def(name = "v_pl850_fc", units = "m s^-1",
                      dim = list(stationdim,memberdim,timedim),
                      missval = fillvalue, longname = dlname,
                      prec="single")

dlname <- "interpolated specific humidity at pressure level 850 hPa ensemble forecast"
q_fc_def <- ncvar_def(name = "q_pl850_fc", units = "kg kg^-1",
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
ncfile_name <- paste0(data_dir,"data_aux_pl850_interpolated_00UTC.nc")
ncout <- nc_create(ncfile_name,
                   list(u_fc_def, v_fc_def, q_fc_def, alt_def, lat_def, lon_def, id_def, location_def),
                   force_v4=T)

# put variables
ncvar_put(nc = ncout, varid = u_fc_def, vals = u_pl850_fc_00)
ncvar_put(nc = ncout, varid = v_fc_def, vals = v_pl850_fc_00)
ncvar_put(nc = ncout, varid = q_fc_def, vals = q_pl850_fc_00)
ncvar_put(nc = ncout, varid = alt_def, vals = st_meta$ALTITUDE)
ncvar_put(nc = ncout, varid = lat_def, vals = st_meta$LATITUDE)
ncvar_put(nc = ncout, varid = lon_def, vals = st_meta$LONGITUDE)
ncvar_put(nc = ncout, varid = id_def, vals = st_meta$STATIONS_ID)
ncvar_put(nc = ncout, varid = location_def, vals = st_meta$LOCATION)

nc_close(ncout)


## 12 UTC file

# define dimensions
rm(timedim)
timedim <- ncdim_def(name = "time", vals = as.integer(time_12),
                     units = "seconds since 1970-01-01 00:00 UTC",
                     longname = "valid time of forecasts and observations, UTC")


# define variables
dlname <- "interpolated u wind component at pressure level 850 hPa ensemble forecast"
u_fc_def <- ncvar_def(name = "u_pl850_fc", units = "m s^-1",
                      dim = list(stationdim,memberdim,timedim),
                      missval = fillvalue, longname = dlname,
                      prec="single")

dlname <- "interpolated v wind component at pressure level 850 hPa ensemble forecast"
v_fc_def <- ncvar_def(name = "v_pl850_fc", units = "m s^-1",
                      dim = list(stationdim,memberdim,timedim),
                      missval = fillvalue, longname = dlname,
                      prec="single")

dlname <- "interpolated specific humidity at pressure level 850 hPa ensemble forecast"
q_fc_def <- ncvar_def(name = "q_pl850_fc", units = "kg kg^-1",
                      dim = list(stationdim,memberdim,timedim),
                      missval = fillvalue, longname = dlname,
                      prec="single")


## create nc file
ncfile_name <- paste0(data_dir,"data_aux_pl850_interpolated_12UTC.nc")
ncout <- nc_create(ncfile_name,
                   list(u_fc_def, v_fc_def, q_fc_def, alt_def, lat_def, lon_def, id_def, location_def),
                   force_v4=T)

# put variables
ncvar_put(nc = ncout, varid = u_fc_def, vals = u_pl850_fc_12)
ncvar_put(nc = ncout, varid = v_fc_def, vals = v_pl850_fc_12)
ncvar_put(nc = ncout, varid = q_fc_def, vals = q_pl850_fc_12)
ncvar_put(nc = ncout, varid = alt_def, vals = st_meta$ALTITUDE)
ncvar_put(nc = ncout, varid = lat_def, vals = st_meta$LATITUDE)
ncvar_put(nc = ncout, varid = lon_def, vals = st_meta$LONGITUDE)
ncvar_put(nc = ncout, varid = id_def, vals = st_meta$STATIONS_ID)
ncvar_put(nc = ncout, varid = location_def, vals = st_meta$LOCATION)

nc_close(ncout)

## -------------------------------------- ###

### surface data ###

nc <- nc_open(paste0(data_dir,"data_aux_surface_interpolated.nc"))
nc

# cape_fc, sp_fc, tcc_fc

cape_fc <- ncvar_get(nc, "cape_fc")
sp_fc <- ncvar_get(nc, "sp_fc")
tcc_fc <- ncvar_get(nc, "tcc_fc")

station_alt <- ncvar_get(nc, "station_alt")
station_lat <- ncvar_get(nc, "station_lat")
station_lon <- ncvar_get(nc, "station_lon")
station_id <- ncvar_get(nc, "station_id")
station_loc <- ncvar_get(nc, "station_loc")

time <- ncvar_get(nc, "time")
station <- ncvar_get(nc, "station")

nc_close(nc)

validtime <- as.POSIXct(time, tz = "UTC", origin = "1970-01-01 00:00")
validtime_hour <- format(validtime, "%H")

ind00 <- which(validtime_hour == "00")
ind12 <- which(validtime_hour == "12")

cape_fc_00 <- cape_fc[,,ind00]
cape_fc_12 <- cape_fc[,,ind12]

sp_fc_00 <- sp_fc[,,ind00]
sp_fc_12 <- sp_fc[,,ind12]

tcc_fc_00 <- tcc_fc[,,ind00]
tcc_fc_12 <- tcc_fc[,,ind12]

time_00 <- time[ind00]
time_12 <- time[ind12]

rm(time, cape_fc, sp_fc, tcc_fc)

## write to nc files

## 00 UTC file

# define dimensions
stationdim <- ncdim_def("station", "station_ID", station_id)
memberdim <- ncdim_def("member", "member_number", as.integer(1:50))
timedim <- ncdim_def(name = "time", vals = as.integer(time_00),
                     units = "seconds since 1970-01-01 00:00 UTC",
                     longname = "valid time of forecasts and observations, UTC")

# define variables
fillvalue <- NA

dlname <- "interpolated cape ensemble forecast"
cape_fc_def <- ncvar_def(name = "cape_fc", units = "J kg^-1",
                         dim = list(stationdim,memberdim,timedim),
                         missval = fillvalue, longname = dlname,
                         prec="single")

dlname <- "interpolated surface pressure ensemble forecast"
sp_fc_def <- ncvar_def(name = "sp_fc", units = "Pa",
                       dim = list(stationdim,memberdim,timedim),
                       missval = fillvalue, longname = dlname,
                       prec="single")

dlname <- "interpolated total cloud cover ensemble forecast"
tcc_fc_def <- ncvar_def(name = "tcc_fc", units = "%",
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
ncfile_name <- paste0(data_dir,"data_aux_surface_interpolated_00UTC.nc")
ncout <- nc_create(ncfile_name,
                   list(cape_fc_def, sp_fc_def, tcc_fc_def, alt_def, lat_def, lon_def, id_def, location_def),
                   force_v4=T)

# put variables
ncvar_put(nc = ncout, varid = cape_fc_def, vals = cape_fc_00)
ncvar_put(nc = ncout, varid = sp_fc_def, vals = sp_fc_00)
ncvar_put(nc = ncout, varid = tcc_fc_def, vals = tcc_fc_00)
ncvar_put(nc = ncout, varid = alt_def, vals = st_meta$ALTITUDE)
ncvar_put(nc = ncout, varid = lat_def, vals = st_meta$LATITUDE)
ncvar_put(nc = ncout, varid = lon_def, vals = st_meta$LONGITUDE)
ncvar_put(nc = ncout, varid = id_def, vals = st_meta$STATIONS_ID)
ncvar_put(nc = ncout, varid = location_def, vals = st_meta$LOCATION)

nc_close(ncout)


## 12 UTC file

# define dimensions
rm(timedim)
timedim <- ncdim_def(name = "time", vals = as.integer(time_12),
                     units = "seconds since 1970-01-01 00:00 UTC",
                     longname = "valid time of forecasts and observations, UTC")


# define variables
dlname <- "interpolated cape ensemble forecast"
cape_fc_def <- ncvar_def(name = "cape_fc", units = "J kg^-1",
                         dim = list(stationdim,memberdim,timedim),
                         missval = fillvalue, longname = dlname,
                         prec="single")

dlname <- "interpolated surface pressure ensemble forecast"
sp_fc_def <- ncvar_def(name = "sp_fc", units = "Pa",
                       dim = list(stationdim,memberdim,timedim),
                       missval = fillvalue, longname = dlname,
                       prec="single")

dlname <- "interpolated total cloud cover ensemble forecast"
tcc_fc_def <- ncvar_def(name = "tcc_fc", units = "%",
                        dim = list(stationdim,memberdim,timedim),
                        missval = fillvalue, longname = dlname,
                        prec="single")


## create nc file
ncfile_name <- paste0(data_dir,"data_aux_surface_interpolated_12UTC.nc")
ncout <- nc_create(ncfile_name,
                   list(cape_fc_def, sp_fc_def, tcc_fc_def, alt_def, lat_def, lon_def, id_def, location_def),
                   force_v4=T)

# put variables
ncvar_put(nc = ncout, varid = cape_fc_def, vals = cape_fc_12)
ncvar_put(nc = ncout, varid = sp_fc_def, vals = sp_fc_12)
ncvar_put(nc = ncout, varid = tcc_fc_def, vals = tcc_fc_12)
ncvar_put(nc = ncout, varid = alt_def, vals = st_meta$ALTITUDE)
ncvar_put(nc = ncout, varid = lat_def, vals = st_meta$LATITUDE)
ncvar_put(nc = ncout, varid = lon_def, vals = st_meta$LONGITUDE)
ncvar_put(nc = ncout, varid = id_def, vals = st_meta$STATIONS_ID)
ncvar_put(nc = ncout, varid = location_def, vals = st_meta$LOCATION)

nc_close(ncout)
