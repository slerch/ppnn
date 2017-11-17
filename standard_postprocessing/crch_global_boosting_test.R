## simpler: use data from 2015 to estimate model for 2016 (fixed training period)

## note: normalization not required, happens automatically in crch.boost?

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

library(scoringRules)
library(crch)

## add some more data 

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




date1 <- as.POSIXct("2016-01-01 00:00", tz = "UTC", origin = "1970-01-01 00:00")

# extract start and end date of rolling training period
enddate <- date1 - 2*24*3600
startdate <- enddate - (365-1)*24*3600

vhour <- "00"
dates_validhour <- format(dates, '%H')
ind_training <- which(dates_validhour == vhour & dates >= startdate & dates <= enddate)

# extract forecast data
ensfc_train <- fcdata[,,ind_training]
ensfc_train_mean <- c(apply(ensfc_train, c(1,3), mean))
ensfc_train_var <- c(apply(ensfc_train, c(1,3), var))
cape_mean_train <- c(cape_mean[,ind_training])
cape_var_train <- c(cape_var[,ind_training])
sp_mean_train <- c(sp_mean[,ind_training])
sp_var_train <- c(sp_var[,ind_training])
tcc_mean_train <- c(tcc_mean[,ind_training])
tcc_var_train <- c(tcc_var[,ind_training])
u_pl850_mean_train <- c(u_pl850_mean[,ind_training])
u_pl850_var_train <- c(u_pl850_var[,ind_training])
v_pl850_mean_train <- c(v_pl850_mean[,ind_training])
v_pl850_var_train <- c(v_pl850_var[,ind_training])
q_pl850_mean_train <- c(q_pl850_mean[,ind_training])
q_pl850_var_train <- c(q_pl850_var[,ind_training])
u_pl500_mean_train <- c(u_pl500_mean[,ind_training])
u_pl500_var_train <- c(u_pl500_var[,ind_training])
v_pl500_mean_train <- c(v_pl500_mean[,ind_training])
v_pl500_var_train <- c(v_pl500_var[,ind_training])
gh_pl500_mean_train <- c(gh_pl500_mean[,ind_training])
gh_pl500_var_train <- c(gh_pl500_var[,ind_training])
sshf_mean_train <- c(sshf_mean[,ind_training])
sshf_var_train <- c(sshf_var[,ind_training])
slhf_mean_train <- c(slhf_mean[,ind_training])
slhf_var_train <- c(slhf_var[,ind_training])
u10_mean_train <- c(u10_mean[,ind_training])
u10_var_train <- c(u10_var[,ind_training])
v10_mean_train <- c(v10_mean[,ind_training])
v10_var_train <- c(v10_var[,ind_training])
ssr_mean_train <- c(ssr_mean[,ind_training])
ssr_var_train <- c(ssr_var[,ind_training])
str_mean_train <- c(str_mean[,ind_training])
str_var_train <- c(str_var[,ind_training])
d2m_mean_train <- c(d2m_mean[,ind_training])
d2m_var_train <- c(d2m_var[,ind_training])
sm_mean_train <- c(sm_mean[,ind_training])
sm_var_train <- c(sm_var[,ind_training])


# extract corresponding observation data
obs_train <- c(obsdata[,ind_training])

if(anyNA(obs_train)){
  ind_notNA <- which(!is.na(obs_train))
  obs_train <- obs_train[ind_notNA]
  ensfc_train_mean <- ensfc_train_mean[ind_notNA]
  ensfc_train_var <- ensfc_train_var[ind_notNA]
  cape_mean_train <- cape_mean_train[ind_notNA]
  cape_var_train <- cape_var_train[ind_notNA]
  sp_mean_train <- sp_mean_train[ind_notNA]
  sp_var_train <- sp_var_train[ind_notNA]
  tcc_mean_train <- tcc_mean_train[ind_notNA]
  tcc_var_train <- tcc_var_train[ind_notNA]
  u_pl850_mean_train <- u_pl850_mean_train[ind_notNA]
  u_pl850_var_train <- u_pl850_var_train[ind_notNA]
  v_pl850_mean_train <- v_pl850_mean_train[ind_notNA]
  v_pl850_var_train <- v_pl850_var_train[ind_notNA]
  q_pl850_mean_train <- q_pl850_mean_train[ind_notNA]
  q_pl850_var_train <- q_pl850_var_train[ind_notNA]
  u_pl500_mean_train <- u_pl500_mean_train[ind_notNA]
  u_pl500_var_train <- u_pl500_var_train[ind_notNA]
  v_pl500_mean_train <- v_pl500_mean_train[ind_notNA]
  v_pl500_var_train <- v_pl500_var_train[ind_notNA]
  gh_pl500_mean_train <- gh_pl500_mean_train[ind_notNA]
  gh_pl500_var_train <- gh_pl500_var_train[ind_notNA]
  sshf_mean_train <- sshf_mean_train[ind_notNA]
  sshf_var_train <- sshf_var_train[ind_notNA]
  slhf_mean_train <- slhf_mean_train[ind_notNA]
  slhf_var_train <- slhf_var_train[ind_notNA]
  u10_mean_train <- u10_mean_train[ind_notNA]
  u10_var_train <- u10_var_train[ind_notNA]
  v10_mean_train <- v10_mean_train[ind_notNA]
  v10_var_train <- v10_var_train[ind_notNA]
  ssr_mean_train <- ssr_mean_train[ind_notNA]
  ssr_var_train <- ssr_var_train[ind_notNA]
  str_mean_train <- str_mean_train[ind_notNA]
  str_var_train <- str_var_train[ind_notNA]
  d2m_mean_train <- d2m_mean_train[ind_notNA]
  d2m_var_train <- d2m_var_train[ind_notNA]
  sm_mean_train <- sm_mean_train[ind_notNA]
  sm_var_train <- sm_var_train[ind_notNA]
}

data_train <- as.data.frame(cbind(obs_train, ensfc_train_mean, ensfc_train_var,
                                  cape_mean_train, cape_var_train,
                                  sp_mean_train, sp_var_train,
                                  tcc_mean_train, tcc_var_train,
                                  u_pl850_mean_train, u_pl850_var_train,
                                  v_pl850_mean_train, v_pl850_var_train,
                                  q_pl850_mean_train, q_pl850_var_train,
                                  u_pl500_mean_train, u_pl500_var_train,
                                  v_pl500_mean_train, v_pl500_var_train,
                                  gh_pl500_mean_train, gh_pl500_var_train,
                                  sshf_mean_train, sshf_var_train,
                                  slhf_mean_train, slhf_var_train,
                                  u10_mean_train, u10_var_train,
                                  v10_mean_train, v10_var_train,
                                  ssr_mean_train, ssr_var_train,
                                  str_mean_train, str_var_train,
                                  d2m_mean_train, d2m_var_train,
                                  sm_mean_train, sm_var_train))
names(data_train) <- c("obs", "ens_mean", "ens_var",
                       "cape_mean", "cape_var", 
                       "sp_mean", "sp_var",
                       "tcc_mean", "tcc_var",
                       "u_pl850_mean", "u_pl850_var",
                       "v_pl850_mean", "v_pl850_var",
                       "q_pl850_mean", "q_pl850_var",
                       "u_pl500_mean", "u_pl500_var",
                       "v_pl500_mean", "v_pl500_var",
                       "gh_pl500_mean", "gh_pl500_var",
                       "sshf_mean", "sshf_var",
                       "slhf_mean", "slhf_var",
                       "u10_mean", "u10_var",
                       "v10_mean", "v10_var",
                       "ssr_mean", "ssr_var",
                       "str_mean", "str_var",
                       "d2m_mean", "d2m_var",
                       "sm_mean", "sm_var")

# model from before
crch_model <- crch(obs ~ ens_mean  | ens_var,
                   data = data_train,
                   dist = "gaussian",
                   link.scale = "log")

# just add some variables
crch_model2 <- crch(obs ~ ens_mean + cape_mean + sp_mean + tcc_mean | ens_var + cape_var + sp_var + tcc_var,
                   data = data_train,
                   dist = "gaussian",
                   link.scale = "log")

# boosting original + surface
crch_model3_boost <- crch(obs ~ .|.,
                          data = data_train[,1:9],
                          dist = "gaussian",
                          link.scale = "log",
                          method = "boosting",
                          mstop = "aic")

# boosting original + surface + pl 850
crch_model4_boost <- crch(obs ~ .|.,
                          data = data_train[,1:15],
                          dist = "gaussian",
                          link.scale = "log",
                          method = "boosting",
                          mstop = "aic")

# boosting original + surface + pl 850 + pl500
crch_model5_boost <- crch(obs ~ .|.,
                          data = data_train[,1:21],
                          dist = "gaussian",
                          link.scale = "log",
                          method = "boosting",
                          mstop = "aic")

# boosting original + surface + pl 850 + pl500 + more surface part 1
crch_model6_boost <- crch(obs ~ .|.,
                          data = data_train[,1:29],
                          dist = "gaussian",
                          link.scale = "log",
                          method = "boosting",
                          mstop = "aic")

# boosting with everything there is 
## actually: ignore soil moisture (contains some NaN values)
crch_model7_boost <- crch(obs ~ .|.,
                          data = data_train[,1:35],
                          dist = "gaussian",
                          link.scale = "log",
                          method = "boosting",
                          mstop = "aic")

## evaluation

startdate_eval <- as.POSIXct("2016-01-01 00:00", tz = "UTC", origin = "1970-01-01 00:00")
enddate_eval <- as.POSIXct("2016-12-31 00:00", tz = "UTC", origin = "1970-01-01 00:00")
ind_eval <- which(dates_validhour == vhour & dates >= startdate_eval & dates <= enddate_eval)
ensfc_eval <- fcdata[,,ind_eval]
ensfc_eval_mean <- c(apply(ensfc_eval, c(1,3), mean))
ensfc_eval_var <- c(apply(ensfc_eval, c(1,3), var))
obs_eval <- c(obsdata[,ind_eval])
cape_mean_eval <- c(cape_mean[,ind_eval])
cape_var_eval <- c(cape_var[,ind_eval])
sp_mean_eval <- c(sp_mean[,ind_eval])
sp_var_eval <- c(sp_var[,ind_eval])
tcc_mean_eval <- c(tcc_mean[,ind_eval])
tcc_var_eval <- c(tcc_var[,ind_eval])
u_pl850_mean_eval <- c(u_pl850_mean[,ind_eval])
u_pl850_var_eval <- c(u_pl850_var[,ind_eval])
v_pl850_mean_eval <- c(v_pl850_mean[,ind_eval])
v_pl850_var_eval <- c(v_pl850_var[,ind_eval])
q_pl850_mean_eval <- c(q_pl850_mean[,ind_eval])
q_pl850_var_eval <- c(q_pl850_var[,ind_eval])
u_pl500_mean_eval <- c(u_pl500_mean[,ind_eval])
u_pl500_var_eval <- c(u_pl500_var[,ind_eval])
v_pl500_mean_eval <- c(v_pl500_mean[,ind_eval])
v_pl500_var_eval <- c(v_pl500_var[,ind_eval])
gh_pl500_mean_eval <- c(gh_pl500_mean[,ind_eval])
gh_pl500_var_eval <- c(gh_pl500_var[,ind_eval])
sshf_mean_eval <- c(sshf_mean[,ind_eval])
sshf_var_eval <- c(sshf_var[,ind_eval])
slhf_mean_eval <- c(slhf_mean[,ind_eval])
slhf_var_eval <- c(slhf_var[,ind_eval])
u10_mean_eval <- c(u10_mean[,ind_eval])
u10_var_eval <- c(u10_var[,ind_eval])
v10_mean_eval <- c(v10_mean[,ind_eval])
v10_var_eval <- c(v10_var[,ind_eval])
ssr_mean_eval <- c(ssr_mean[,ind_eval])
ssr_var_eval <- c(ssr_var[,ind_eval])
str_mean_eval <- c(str_mean[,ind_eval])
str_var_eval <- c(str_var[,ind_eval])
d2m_mean_eval <- c(d2m_mean[,ind_eval])
d2m_var_eval <- c(d2m_var[,ind_eval])
sm_mean_eval <- c(sm_mean[,ind_eval])
sm_var_eval <- c(sm_var[,ind_eval])

if(anyNA(obs_eval)){
  ind_notNA <- which(!is.na(obs_eval))
  obs_eval <- obs_eval[ind_notNA]
  ensfc_eval_mean <- ensfc_eval_mean[ind_notNA]
  ensfc_eval_var <- ensfc_eval_var[ind_notNA]
  cape_mean_eval <- cape_mean_eval[ind_notNA]
  cape_var_eval <- cape_var_eval[ind_notNA]
  sp_mean_eval <- sp_mean_eval[ind_notNA]
  sp_var_eval <- sp_var_eval[ind_notNA]
  tcc_mean_eval <- tcc_mean_eval[ind_notNA]
  tcc_var_eval <- tcc_var_eval[ind_notNA]
  u_pl850_mean_eval <- u_pl850_mean_eval[ind_notNA]
  u_pl850_var_eval <- u_pl850_var_eval[ind_notNA]
  v_pl850_mean_eval <- v_pl850_mean_eval[ind_notNA]
  v_pl850_var_eval <- v_pl850_var_eval[ind_notNA]
  q_pl850_mean_eval <- q_pl850_mean_eval[ind_notNA]
  q_pl850_var_eval <- q_pl850_var_eval[ind_notNA]
  u_pl500_mean_eval <- u_pl500_mean_eval[ind_notNA]
  u_pl500_var_eval <- u_pl500_var_eval[ind_notNA]
  v_pl500_mean_eval <- v_pl500_mean_eval[ind_notNA]
  v_pl500_var_eval <- v_pl500_var_eval[ind_notNA]
  gh_pl500_mean_eval <- gh_pl500_mean_eval[ind_notNA]
  gh_pl500_var_eval <- gh_pl500_var_eval[ind_notNA]
  sshf_mean_eval <- sshf_mean_eval[ind_notNA]
  sshf_var_eval <- sshf_var_eval[ind_notNA]
  slhf_mean_eval <- slhf_mean_eval[ind_notNA]
  slhf_var_eval <- slhf_var_eval[ind_notNA]
  u10_mean_eval <- u10_mean_eval[ind_notNA]
  u10_var_eval <- u10_var_eval[ind_notNA]
  v10_mean_eval <- v10_mean_eval[ind_notNA]
  v10_var_eval <- v10_var_eval[ind_notNA]
  ssr_mean_eval <- ssr_mean_eval[ind_notNA]
  ssr_var_eval <- ssr_var_eval[ind_notNA]
  str_mean_eval <- str_mean_eval[ind_notNA]
  str_var_eval <- str_var_eval[ind_notNA]
  d2m_mean_eval <- d2m_mean_eval[ind_notNA]
  d2m_var_eval <- d2m_var_eval[ind_notNA]
  sm_mean_eval <- sm_mean_eval[ind_notNA]
  sm_var_eval <- sm_var_eval[ind_notNA]
}
data_eval <- as.data.frame(cbind(obs_eval, ensfc_eval_mean, ensfc_eval_var,
                                 cape_mean_eval, cape_var_eval,
                                 sp_mean_eval, sp_var_eval,
                                 tcc_mean_eval, tcc_var_eval,
                                 u_pl850_mean_eval, u_pl850_var_eval,
                                 v_pl850_mean_eval, v_pl850_var_eval,
                                 q_pl850_mean_eval, q_pl850_var_eval,
                                 u_pl500_mean_eval, u_pl500_var_eval,
                                 v_pl500_mean_eval, v_pl500_var_eval,
                                 gh_pl500_mean_eval, gh_pl500_var_eval,
                                 sshf_mean_eval, sshf_var_eval,
                                 slhf_mean_eval, slhf_var_eval,
                                 u10_mean_eval, u10_var_eval,
                                 v10_mean_eval, v10_var_eval,
                                 ssr_mean_eval, ssr_var_eval,
                                 str_mean_eval, str_var_eval,
                                 d2m_mean_eval, d2m_var_eval,
                                 sm_mean_eval, sm_var_eval))
names(data_eval) <- c("obs", "ens_mean", "ens_var",
                      "cape_mean", "cape_var", 
                      "sp_mean", "sp_var",
                      "tcc_mean", "tcc_var",
                      "u_pl850_mean", "u_pl850_var",
                      "v_pl850_mean", "v_pl850_var",
                      "q_pl850_mean", "q_pl850_var",
                      "u_pl500_mean", "u_pl500_var",
                      "v_pl500_mean", "v_pl500_var",
                      "gh_pl500_mean", "gh_pl500_var",
                      "sshf_mean", "sshf_var",
                      "slhf_mean", "slhf_var",
                      "u10_mean", "u10_var",
                      "v10_mean", "v10_var",
                      "ssr_mean", "ssr_var",
                      "str_mean", "str_var",
                      "d2m_mean", "d2m_var",
                      "sm_mean", "sm_var")

loc1 <- as.numeric(predict(crch_model, data_eval, type = "location"))
sc1 <- as.numeric(predict(crch_model, data_eval, type = "scale"))
mean(crps_norm(y = data_eval$obs, mean = loc1, sd = sc1)) # 1.013909

loc2 <- as.numeric(predict(crch_model2, data_eval, type = "location"))
sc2 <- as.numeric(predict(crch_model2, data_eval, type = "scale"))
mean(crps_norm(y = data_eval$obs, mean = loc2, sd = sc2)) # 0.9961115

loc3 <- as.numeric(predict(crch_model3_boost, data_eval, type = "location"))
sc3 <- as.numeric(predict(crch_model3_boost, data_eval, type = "scale"))
mean(crps_norm(y = data_eval$obs, mean = loc3, sd = sc3)) # 0.9925455

loc4 <- as.numeric(predict(crch_model4_boost, data_eval, type = "location"))
sc4 <- as.numeric(predict(crch_model4_boost, data_eval, type = "scale"))
mean(crps_norm(y = data_eval$obs, mean = loc4, sd = sc4)) # 0.9961843 (worse than model 3)

loc5 <- as.numeric(predict(crch_model5_boost, data_eval, type = "location"))
sc5 <- as.numeric(predict(crch_model5_boost, data_eval, type = "scale"))
mean(crps_norm(y = data_eval$obs, mean = loc5, sd = sc5)) # 0.9989608

loc6 <- as.numeric(predict(crch_model6_boost, data_eval, type = "location"))
sc6 <- as.numeric(predict(crch_model6_boost, data_eval, type = "scale"))
mean(crps_norm(y = data_eval$obs, mean = loc6, sd = sc6)) # 0.988455

loc7 <- as.numeric(predict(crch_model7_boost, data_eval, type = "location"))
sc7 <- as.numeric(predict(crch_model7_boost, data_eval, type = "scale"))
mean(crps_norm(y = data_eval$obs, mean = loc7, sd = sc7)) # 0.9765516

## test mstop = "max", "bic", "cv"

## local model implementation: loop over stations, use predict() directly within loop and save out of sample CRPS values