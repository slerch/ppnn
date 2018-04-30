## prepare data set as Rdata file for easier loading and use on cluster 
## (which has no nc4 support...)

rm(list=ls())

## load data
data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
library(ncdf4)
nc <- nc_open(paste0(data_dir, "data_interpolated_00UTC.nc"))
fcdata <- ncvar_get(nc, "t2m_fc") 
obsdata <- ncvar_get(nc, "t2m_obs")
dates <- as.POSIXct(ncvar_get(nc, "time"), origin = "1970-01-01 00:00", tz = "UTC")
stations <- ncvar_get(nc, "station")

station_metadata <- list()
station_metadata$altitudes <- ncvar_get(nc, "station_alt")
station_metadata$latitudes <- ncvar_get(nc, "station_lat")
station_metadata$longitudes <- ncvar_get(nc, "station_lon")
station_metadata$locations <- ncvar_get(nc, "station_loc")
nc_close(nc)

## add covariate data 

## surface data 

nc <- nc_open(paste0(data_dir, "data_aux_surface_interpolated_00UTC.nc"))

cape <- ncvar_get(nc, "cape_fc")
sp <- ncvar_get(nc, "sp_fc") # surface pressure
tcc <- ncvar_get(nc, "tcc_fc")

nc_close(nc)


# keep only mean and var

cape_mean <- apply(cape, c(1,3), mean)
cape_var <- apply(cape, c(1,3), var)
rm(cape)

sp_mean <- apply(sp, c(1,3), mean)
sp_var <- apply(sp, c(1,3), var)
rm(sp)

tcc_mean <- apply(tcc, c(1,3), mean)
tcc_var <- apply(tcc, c(1,3), var)
rm(tcc)


## data from pressure level 850 hPa

nc <- nc_open(paste0(data_dir, "data_aux_pl850_interpolated_00UTC.nc"))

u_pl850 <- ncvar_get(nc, "u_pl850_fc")
v_pl850 <- ncvar_get(nc, "v_pl850_fc")
q_pl850 <- ncvar_get(nc, "q_pl850_fc") # specific humidity

nc_close(nc)

# keep only mean and var

u_pl850_mean <- apply(u_pl850, c(1,3), mean)
u_pl850_var <- apply(u_pl850, c(1,3), var)
rm(u_pl850)

v_pl850_mean <- apply(v_pl850, c(1,3), mean)
v_pl850_var <- apply(v_pl850, c(1,3), var)
rm(v_pl850)

q_pl850_mean <- apply(q_pl850, c(1,3), mean)
q_pl850_var <- apply(q_pl850, c(1,3), var)
rm(q_pl850)


## data from pressure level 500 hPa

nc <- nc_open(paste0(data_dir, "data_aux_pl500_interpolated_00UTC.nc"))

u_pl500 <- ncvar_get(nc, "u_pl500_fc")
v_pl500 <- ncvar_get(nc, "v_pl500_fc")
gh_pl500 <- ncvar_get(nc, "gh_pl500_fc") # geopotential height

nc_close(nc)

# keep only mean and var

u_pl500_mean <- apply(u_pl500, c(1,3), mean)
u_pl500_var <- apply(u_pl500, c(1,3), var)
rm(u_pl500)

v_pl500_mean <- apply(v_pl500, c(1,3), mean)
v_pl500_var <- apply(v_pl500, c(1,3), var)
rm(v_pl500)

gh_pl500_mean <- apply(gh_pl500, c(1,3), mean)
gh_pl500_var <- apply(gh_pl500, c(1,3), var)
rm(gh_pl500)


## more surface variables - part 1 of 2

nc <- nc_open(paste0(data_dir, "data_aux_surface_more_interpolated_part1_00UTC.nc"))

sshf <- ncvar_get(nc, "sshf_fc") # sensible heat flux
slhf <- ncvar_get(nc, "slhf_fc") # latent heat flux
u10 <- ncvar_get(nc, "u10_fc") # u wind 
v10 <- ncvar_get(nc, "v10_fc") # u wind 

nc_close(nc)

# keep only mean and var

sshf_mean <- apply(sshf, c(1,3), mean)
sshf_var <- apply(sshf, c(1,3), var)
rm(sshf)

slhf_mean <- apply(slhf, c(1,3), mean)
slhf_var <- apply(slhf, c(1,3), var)
rm(slhf)

u10_mean <- apply(u10, c(1,3), mean)
u10_var <- apply(u10, c(1,3), var)
rm(u10)

v10_mean <- apply(v10, c(1,3), mean)
v10_var <- apply(v10, c(1,3), var)
rm(v10)


## more surface variables - part 2 of 2

nc <- nc_open(paste0(data_dir, "data_aux_surface_more_interpolated_part2_00UTC.nc"))

ssr <- ncvar_get(nc, "ssr_fc") # net short wave radiation flux
str <- ncvar_get(nc, "str_fc") # net long wave radiation flux 
d2m <- ncvar_get(nc, "d2m_fc") # 2m dew point temperature
sm <- ncvar_get(nc, "sm_fc") # soil moisture

nc_close(nc)

ssr_mean <- apply(ssr, c(1,3), mean)
ssr_var <- apply(ssr, c(1,3), var)
rm(ssr)

str_mean <- apply(str, c(1,3), mean)
str_var <- apply(str, c(1,3), var)
rm(str)

d2m_mean <- apply(d2m, c(1,3), mean)
d2m_var <- apply(d2m, c(1,3), var)
rm(d2m)

sm_mean <- apply(sm, c(1,3), mean)
sm_var <- apply(sm, c(1,3), var)
rm(sm)
## note: sm contains NAs


##
## combine data into a data frame for easier handling 
##

# transform predictors to vector, add date, station and obs

# ens fc of t2m
t2m_mean <- apply(fcdata, c(1,3), mean)
t2m_var <- apply(fcdata, c(1,3), var)

# observations
obs <- c(obsdata)

rm(fcdata, obsdata)

# compute date and station vector
# NOTE: matrices have stations in row and date in column
# c(.) puts together columns into a vector
# here is a trivial example of matrix to vector transformations:
# mm <- matrix(paste0("st",rep(1:3, each = 2)), nrow = 3, byrow = T)
# mm[,1] <- paste0("day1-", mm[,1])
# mm[,2] <- paste0("day2-", mm[,2])
# mm
# c(mm)
# rep(c("day1", "day2"), each = 3)
# rep(paste0("st", 1:3), 2)
# rm(mm)

dates_vec <- rep(dates, each = 537)
stations_vec <- rep(stations, length(dates))

data <- data.frame(dates_vec,
                   stations_vec,
                   obs, 
                   c(t2m_mean), 
                   c(t2m_var),
                   c(cape_mean), 
                   c(cape_var),
                   c(sp_mean), 
                   c(sp_var),
                   c(tcc_mean),
                   c(tcc_var),
                   c(u_pl850_mean), 
                   c(u_pl850_var),
                   c(v_pl850_mean),
                   c(v_pl850_var),
                   c(q_pl850_mean), 
                   c(q_pl850_var),
                   c(u_pl500_mean), 
                   c(u_pl500_var),
                   c(v_pl500_mean), 
                   c(v_pl500_var),
                   c(gh_pl500_mean), 
                   c(gh_pl500_var),
                   c(sshf_mean), 
                   c(sshf_var),
                   c(slhf_mean), 
                   c(slhf_var),
                   c(u10_mean), 
                   c(u10_var),
                   c(v10_mean), 
                   c(v10_var),
                   c(ssr_mean), 
                   c(ssr_var),
                   c(str_mean), 
                   c(str_var),
                   c(d2m_mean), 
                   c(d2m_var),
                   c(sm_mean), 
                   c(sm_var)
)


# modify names of columns
names(data)[4:ncol(data)] <- unlist(lapply(names(data)[4:ncol(data)],
                                          function(z) substr(z, 3, nchar(z) -1 )))

names(data)[1:2] <- c("date", "station")

# check output
head(data)

# save output
save(data,
     file = paste0(data_dir, "data_all.Rdata"))
