rm(list=ls())
library(ncdf4)

data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
nc <- nc_open(paste0(data_dir,"data_interpolated_corrected.nc"))
nc

t2m_fc <- ncvar_get(nc, "t2m_fc")
# dim(t2m_fc)
# hist(t2m_fc)

t2m_obs <- ncvar_get(nc, "t2m_obs")
station_alt <- ncvar_get(nc, "station_alt")
station_lat <- ncvar_get(nc, "station_lat")
station_lon <- ncvar_get(nc, "station_lon")
station_id <- ncvar_get(nc, "station_id")
station_loc <- ncvar_get(nc, "station_loc")

time <- ncvar_get(nc, "time")
station <- ncvar_get(nc, "station")

validtime <- as.POSIXct(time, tz = "UTC", origin = "1970-01-01 00:00")
validtime_hour <- format(validtime, "%H")
  
ind00 <- which(validtime_hour == "00")
ind12 <- which(validtime_hour == "12")

# hist(t2m_fc[,,ind12])

t2m_fc_00 <- t2m_fc[,,ind00]
t2m_fc_12 <- t2m_fc[,,ind12]
t2m_obs_00 <- t2m_obs[,ind00]
t2m_obs_12 <- t2m_obs[,ind12]
time_00 <- time[ind00]
time_12 <- time[ind12]

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


## 00 UTC file

## create nc file
ncfile_name <- paste0(data_dir,"data_interpolated_00UTC.nc")
ncout <- nc_create(ncfile_name,
                   list(t2m_fc_def, t2m_obs_def, alt_def, lat_def, lon_def, id_def, location_def),
                   force_v4=T)

# put variables
ncvar_put(nc = ncout, varid = t2m_fc_def, vals = t2m_fc_00)
ncvar_put(nc = ncout, varid = t2m_obs_def, vals = t2m_obs_00)
ncvar_put(nc = ncout, varid = alt_def, vals = station_alt)
ncvar_put(nc = ncout, varid = lat_def, vals = station_lat)
ncvar_put(nc = ncout, varid = lon_def, vals = station_lon)
ncvar_put(nc = ncout, varid = id_def, vals = station_id)
ncvar_put(nc = ncout, varid = location_def, vals = station_loc)

nc_close(ncout)


## 12 UTC file

# define dimensions
rm(timedim)
timedim <- ncdim_def(name = "time", vals = as.integer(time_12),
                     units = "seconds since 1970-01-01 00:00 UTC",
                     longname = "valid time of forecasts and observations, UTC")


# define variables
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


## create nc file
ncfile_name <- paste0(data_dir,"data_interpolated_12UTC.nc")
ncout <- nc_create(ncfile_name,
                   list(t2m_fc_def, t2m_obs_def, alt_def, lat_def, lon_def, id_def, location_def),
                   force_v4=T)

# put variables
ncvar_put(nc = ncout, varid = t2m_fc_def, vals = t2m_fc_12)
ncvar_put(nc = ncout, varid = t2m_obs_def, vals = t2m_obs_12)
ncvar_put(nc = ncout, varid = alt_def, vals = station_alt)
ncvar_put(nc = ncout, varid = lat_def, vals = station_lat)
ncvar_put(nc = ncout, varid = lon_def, vals = station_lon)
ncvar_put(nc = ncout, varid = id_def, vals = station_id)
ncvar_put(nc = ncout, varid = location_def, vals = station_loc)

nc_close(ncout)
