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

## part 1 of aux_surface_more

nc <- nc_open(paste0(data_dir,"data_aux_surface_more_interpolated_part1.nc"))
nc

# sshf_fc, slhf_fc, u10_fc, v10_fc

sshf_fc <- ncvar_get(nc, "sshf_fc")
slhf_fc <- ncvar_get(nc, "slhf_fc")
u10_fc <- ncvar_get(nc, "u10_fc")
v10_fc <- ncvar_get(nc, "v10_fc")

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

sshf_fc_00 <- sshf_fc[,,ind00]
sshf_fc_12 <- sshf_fc[,,ind12]

slhf_fc_00 <- slhf_fc[,,ind00]
slhf_fc_12 <- slhf_fc[,,ind12]

u10_fc_00 <- u10_fc[,,ind00]
u10_fc_12 <- u10_fc[,,ind12]

v10_fc_00 <- v10_fc[,,ind00]
v10_fc_12 <- v10_fc[,,ind12]

time_00 <- time[ind00]
time_12 <- time[ind12]

rm(time, sshf_fc, slhf_fc, u10_fc, v10_fc)


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

dlname <- "interpolated sensible heat flux ensemble forecast"
sshf_fc_def <- ncvar_def(name = "sshf_fc", units = "W m^-2",
                         dim = list(stationdim,memberdim,timedim),
                         missval = fillvalue, longname = dlname,
                         prec="single")

dlname <- "interpolated latent heat flux ensemble forecast"
slhf_fc_def <- ncvar_def(name = "slhf_fc", units = "W m^-2",
                         dim = list(stationdim,memberdim,timedim),
                         missval = fillvalue, longname = dlname,
                         prec="single")

dlname <- "interpolated u-wind ensemble forecast"
u10_fc_def <- ncvar_def(name = "u10_fc", units = "m s^-1",
                        dim = list(stationdim,memberdim,timedim),
                        missval = fillvalue, longname = dlname,
                        prec="single")

dlname <- "interpolated v-wind ensemble forecast"
v10_fc_def <- ncvar_def(name = "v10_fc", units = "m s^-1",
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
ncfile_name <- paste0(data_dir,"data_aux_surface_more_interpolated_part1_00UTC.nc")
ncout <- nc_create(ncfile_name,
                   list(sshf_fc_def, slhf_fc_def, u10_fc_def, v10_fc_def, alt_def, lat_def, lon_def, id_def, location_def),
                   force_v4=T)

# put variables
ncvar_put(nc = ncout, varid = sshf_fc_def, vals = sshf_fc_00)
ncvar_put(nc = ncout, varid = slhf_fc_def, vals = slhf_fc_00)
ncvar_put(nc = ncout, varid = u10_fc_def, vals = u10_fc_00)
ncvar_put(nc = ncout, varid = v10_fc_def, vals = v10_fc_00)
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
fillvalue <- NA

dlname <- "interpolated sensible heat flux ensemble forecast"
sshf_fc_def <- ncvar_def(name = "sshf_fc", units = "W m^-2",
                         dim = list(stationdim,memberdim,timedim),
                         missval = fillvalue, longname = dlname,
                         prec="single")

dlname <- "interpolated latent heat flux ensemble forecast"
slhf_fc_def <- ncvar_def(name = "slhf_fc", units = "W m^-2",
                         dim = list(stationdim,memberdim,timedim),
                         missval = fillvalue, longname = dlname,
                         prec="single")

dlname <- "interpolated u-wind ensemble forecast"
u10_fc_def <- ncvar_def(name = "u10_fc", units = "m s^-1",
                        dim = list(stationdim,memberdim,timedim),
                        missval = fillvalue, longname = dlname,
                        prec="single")

dlname <- "interpolated v-wind ensemble forecast"
v10_fc_def <- ncvar_def(name = "v10_fc", units = "m s^-1",
                        dim = list(stationdim,memberdim,timedim),
                        missval = fillvalue, longname = dlname,
                        prec="single")

## create nc file
ncfile_name <- paste0(data_dir,"data_aux_surface_more_interpolated_part1_12UTC.nc")
ncout <- nc_create(ncfile_name,
                   list(sshf_fc_def, slhf_fc_def, u10_fc_def, v10_fc_def, alt_def, lat_def, lon_def, id_def, location_def),
                   force_v4=T)

# put variables
ncvar_put(nc = ncout, varid = sshf_fc_def, vals = sshf_fc_12)
ncvar_put(nc = ncout, varid = slhf_fc_def, vals = slhf_fc_12)
ncvar_put(nc = ncout, varid = u10_fc_def, vals = u10_fc_12)
ncvar_put(nc = ncout, varid = v10_fc_def, vals = v10_fc_12)
ncvar_put(nc = ncout, varid = alt_def, vals = st_meta$ALTITUDE)
ncvar_put(nc = ncout, varid = lat_def, vals = st_meta$LATITUDE)
ncvar_put(nc = ncout, varid = lon_def, vals = st_meta$LONGITUDE)
ncvar_put(nc = ncout, varid = id_def, vals = st_meta$STATIONS_ID)
ncvar_put(nc = ncout, varid = location_def, vals = st_meta$LOCATION)

nc_close(ncout)

##
##
##

## part 2 data

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

## part 2 of aux_surface_more

nc <- nc_open(paste0(data_dir,"data_aux_surface_more_interpolated_part2.nc"))
nc


# ssr_fc, str_fc, d2m_fc, sm_fc

ssr_fc <- ncvar_get(nc, "ssr_fc")
str_fc <- ncvar_get(nc, "str_fc")
d2m_fc <- ncvar_get(nc, "d2m_fc")
sm_fc <- ncvar_get(nc, "sm_fc")

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

ssr_fc_00 <- ssr_fc[,,ind00]
ssr_fc_12 <- ssr_fc[,,ind12]

str_fc_00 <- str_fc[,,ind00]
str_fc_12 <- str_fc[,,ind12]

d2m_fc_00 <- d2m_fc[,,ind00]
d2m_fc_12 <- d2m_fc[,,ind12]

sm_fc_00 <- sm_fc[,,ind00]
sm_fc_12 <- sm_fc[,,ind12]

time_00 <- time[ind00]
time_12 <- time[ind12]

rm(time, ssr_fc, str_fc, d2m_fc, sm_fc)

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
ncfile_name <- paste0(data_dir,"data_aux_surface_more_interpolated_part2_00UTC.nc")
ncout <- nc_create(ncfile_name,
                   list(ssr_fc_def, str_fc_def, d2m_fc_def, sm_fc_def, alt_def, lat_def, lon_def, id_def, location_def),
                   force_v4=T)

# ssr_fc, str_fc, d2m_fc, sm_fc

# put variables
ncvar_put(nc = ncout, varid = ssr_fc_def, vals = ssr_fc_00)
ncvar_put(nc = ncout, varid = str_fc_def, vals = str_fc_00)
ncvar_put(nc = ncout, varid = d2m_fc_def, vals = d2m_fc_00)
ncvar_put(nc = ncout, varid = sm_fc_def, vals = sm_fc_00)
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

## create nc file
ncfile_name <- paste0(data_dir,"data_aux_surface_more_interpolated_part2_12UTC.nc")
ncout <- nc_create(ncfile_name,
                   list(ssr_fc_def, str_fc_def, d2m_fc_def, sm_fc_def, alt_def, lat_def, lon_def, id_def, location_def),
                   force_v4=T)

# ssr_fc, str_fc, d2m_fc, sm_fc

# put variables
ncvar_put(nc = ncout, varid = ssr_fc_def, vals = ssr_fc_12)
ncvar_put(nc = ncout, varid = str_fc_def, vals = str_fc_12)
ncvar_put(nc = ncout, varid = d2m_fc_def, vals = d2m_fc_12)
ncvar_put(nc = ncout, varid = sm_fc_def, vals = sm_fc_12)
ncvar_put(nc = ncout, varid = alt_def, vals = st_meta$ALTITUDE)
ncvar_put(nc = ncout, varid = lat_def, vals = st_meta$LATITUDE)
ncvar_put(nc = ncout, varid = lon_def, vals = st_meta$LONGITUDE)
ncvar_put(nc = ncout, varid = id_def, vals = st_meta$STATIONS_ID)
ncvar_put(nc = ncout, varid = location_def, vals = st_meta$LOCATION)

nc_close(ncout)