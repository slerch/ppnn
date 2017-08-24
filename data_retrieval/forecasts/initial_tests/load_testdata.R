## to process in R, downloaded grib files can be converted to netCDF via
##  grib_to_netcdf -o ecmwf_t2m_2016.nc ecmwf_t2m_2016.grib 
## using the grib_to_netcdf function from ecCodes, 
## see https://software.ecmwf.int/wiki/display/ECC/ecCodes+Home

rm(list=ls())

library(ncdf4)

fcdat <- nc_open("/media/sebastian/Elements/Postproc_NN/data/forecasts/ecmwf_t2m_2016.nc")

t2m <- ncvar_get(fcdat, "t2m") 
# t2m[longitude,latitude,number,time] 

time <- ncvar_get(fcdat, "time")
ncatt_get(fcdat, "time")
validtime <- as.POSIXlt(time*3600, origin = "1900-01-01", tz = "UTC")

t2m[1,1,,1]
hist(c(t2m))