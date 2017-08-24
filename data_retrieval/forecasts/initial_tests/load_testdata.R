## to process in R, downloaded grib files can be converted to netCDF via
##  grib_to_netcdf -o ecmwf_t2m_2016.nc ecmwf_t2m_2016.grib 
## using the grib_to_netcdf function from ecCodes, 
## see https://software.ecmwf.int/wiki/display/ECC/ecCodes+Home

rm(list=ls())

library(ncdf4)

fcdat <- nc_open("/media/sebastian/Elements/Postproc_NN/data/forecasts/ecmwf_t2m_2007.nc")

t2m <- ncvar_get(fcdat, "t2m") 
# t2m[longitude,latitude,number,time] 

time <- ncvar_get(fcdat, "time")
ncatt_get(fcdat, "time")
validtime <- as.POSIXlt(time*3600, origin = "1900-01-01", tz = "UTC")

t2m[1,1,,500]
hist(c(t2m))

lon <- ncvar_get(fcdat, "longitude")
lat <- ncvar_get(fcdat, "latitude")

# mean 12UTC temperature forecast for grid point near Heidelberg
t2m_means <- apply(t2m[which(lon == 8.5),which(lat == 49.5),,], 2, mean)
plot(validtime[c(TRUE,FALSE)], t2m_means[c(TRUE,FALSE)] - 273.15, type = "l")
# add 00 UTC mean forecasts
lines(validtime[c(FALSE,TRUE)], t2m_means[c(FALSE,TRUE)] - 273.15, lty = 2)