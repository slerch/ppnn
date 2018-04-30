## extract ensemble forecasts for evaluation purposes

rm(list=ls())

## load data
data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
library(ncdf4)
nc <- nc_open(paste0(data_dir, "data_interpolated_00UTC.nc"))
fcdata <- ncvar_get(nc, "t2m_fc") 
obsdata <- ncvar_get(nc, "t2m_obs")
dates <- as.POSIXct(ncvar_get(nc, "time"), origin = "1970-01-01 00:00", tz = "UTC")
stations <- ncvar_get(nc, "station")

nc_close(nc)

dim(fcdata)

obs <- c(obsdata)
dates_vec <- rep(dates, each = 537)
stations_vec <- rep(stations, length(dates))

library(abind)
t2m <- adrop(fcdata[1,,,drop=FALSE], drop = 1)
dim(t2m)

t2m <- array(dim = c(537*3653,50))
for(dd in 1:length(dates)){
  if(dd %% 100 == 0){print(dd)}
  ind <- which(dates_vec == dates[dd])
  fc_today <- fcdata[,,dd]
  t2m[ind,] <- fc_today
}

dim(t2m)

# # compare to mean computed as in data_preparation.R
# a <- apply(t2m, 1, mean)
# b.matrix <- apply(fcdata, c(1,3), mean)
# b <- c(b.matrix)
# any(a!=b)

data_ens <- data.frame(date = dates_vec,
                   station = stations_vec,
                   obs, 
                   t2m
)

save(data_ens,
     file = paste0(data_dir, "data_ensfc.Rdata"))

