## implement post-processing with crch model in order to implement boosting for EMOS models later on

## note that CRAN version of crch uses ML estimation (not minimum CRPS)

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

# vector of validhours (00/12 UTC) to be able to select only forecast cases 
# with that specific valid hour for training
dates_validhour <- format(dates, '%H')

postproc_crch <- function(vdate, train_length){
  # determine training set
  vhour <- format(vdate, "%H")
  
  # extract start and end date of rolling training period
  enddate <- vdate - 2*24*3600
  startdate <- enddate - (train_length-1)*24*3600
  if(startdate < dates[1]){
    stop("start date of training period before first available date, choose shorter train_length")
  }
  
  # generate index vector for training set, indicating those forecast cases which
  #   - have to same valid hour as vdate (00/12 UTC)
  #   - are in the training periof between startdate and enddate
  ind_training <- which(dates_validhour == vhour & dates >= startdate & dates <= enddate)
  
  # extract forecast data
  ensfc_train <- fcdata[,,ind_training]
  ensfc_train_mean <- c(apply(ensfc_train, c(1,3), mean))
  ensfc_train_var <- c(apply(ensfc_train, c(1,3), var))
  
  # extract corresponding observation data
  obs_train <- c(obsdata[,ind_training])
  
  # find and omit NA values in observations (there are no NAs in forecasts)
  if(anyNA(obs_train)){
    ind_notNA <- which(!is.na(obs_train))
    obs_train <- obs_train[ind_notNA]
    ensfc_train_mean <- ensfc_train_mean[ind_notNA]
    ensfc_train_var <- ensfc_train_var[ind_notNA]
  }
  
  # # determine optimal EMOS coefficients a,b,c,d using minimum CRPS estimation
  # optim_out <- optim(par = c(1,1,1,1), fn = objective_fun,
  #                    gr = gradfun_wrapper,
  #                    ensmean = ensfc_train_mean, ensvar = ensfc_train_var,
  #                    obs = obs_train,
  #                    method = "BFGS")
  
  ## now replaced by crch model estimation
  crch_model <- crch(obs_train ~ ensfc_train_mean | ensfc_train_var,
                     dist = "gaussian",
                     link.scale = "identity")
  
  crch_par <- as.numeric(c(crch_model$coefficients$location,
                           crch_model$coefficients$scale))
  
  return(crch_par)
}


# example for one date
# tt <- as.POSIXct("2013-11-25 00:00", tz = "UTC", origin = "1970-01-01 00:00")
# m <- 25
# par_out <- postproc_crch(vdate = tt, train_length = m)
# 
# ind_today <- which(dates == tt)
# ensfc_today <- fcdata[,,ind_today]
# ensfc_today_mean <- apply(ensfc_today, c(1), mean)
# ensfc_today_var <- apply(ensfc_today, c(1), var)
# 
# loc <- c(cbind(1, ensfc_today_mean) %*% par_out[1:2])
# sc <- c(cbind(1, ensfc_today_var) %*% par_out[3:4]) # note changes compared to my implementation: crch models scale rather than variance
# 
# obs_today <- obsdata[,ind_today]
# 
# if(anyNA(obs_today)){
#   ind_notNA_today <- which(!is.na(obs_today))
#   obs_today <- obs_today[ind_notNA_today]
#   loc <- loc[ind_notNA_today]
#   sc <- sc[ind_notNA_today]
#   ensfc_today <- ensfc_today[ind_notNA_today,]
# }
# 
# crps_n <- crps_norm(obs_today, loc, sc)
# crps_ens <- crps_sample(obs_today, ensfc_today)
# 
# summary(crps_n)
# summary(crps_ens)
# summary(crps_ens - crps_n)
# hist(crps_ens - crps_n)


## estimate model for 2016
# training length: 25 days
m <- 25
# dates for which EMOS models should be fit (here: only for 00 UTC valid times; all dates in 2008-2016)
dates_fit <- seq(as.POSIXct("2016-01-01 00:00", tz = "UTC", origin = "1970-01-01 00:00"), 
                 as.POSIXct("2016-12-31 00:00", tz = "UTC", origin = "1970-01-01 00:00"), 
                 by="+1 day")
# objects to save mean CRPS values to
crpsout_n <- numeric(length(dates_fit))
crpsout_ens <- numeric(length(dates_fit))

for(dd in 1:length(dates_fit)){
  # choose current date
  day <- dates_fit[dd]
  # progress indicator
  if(format(dates_fit[dd], "%d") == "01"){
    format(dates_fit[dd], "%B %Y")
    cat("Starting at", paste(Sys.time()), ":", format(dates_fit[dd], "%B %Y"), "\n"); flush(stdout())
  }
  # apply EMOS function to obtain optimal parameters
  par_out <- postproc_crch(vdate = day, train_length = m)
  
  # extract ensemble forecasts and summary statistics for the day of interest ("day")
  ind_today <- which(dates == day)
  ensfc_today <- fcdata[,,ind_today]
  ensfc_today_mean <- apply(ensfc_today, c(1), mean)
  ensfc_today_var <- apply(ensfc_today, c(1), var)
  
  # use summary statistics and EMOS coefficients to determine mu and sigma
  # (= location and scale)
  # note that the link functions are the same as employed in the objective function
  loc <- c(cbind(1, ensfc_today_mean) %*% par_out[1:2])
  sc <- c(cbind(1, ensfc_today_var) %*% par_out[3:4])
  
  # to evaluate forecasts, find corresponding observations
  obs_today <- obsdata[,ind_today]
  
  # loc, sc and obs_today are vectors, with values for all stations at date "day" 
  # there are some NA values in observations: omit those forecast cases 
  if(anyNA(obs_today)){
    ind_notNA_today <- which(!is.na(obs_today))
    obs_today <- obs_today[ind_notNA_today]
    loc <- loc[ind_notNA_today]
    sc <- sc[ind_notNA_today]
    ensfc_today <- ensfc_today[ind_notNA_today,]
  }
  
  # compute the mean CRPS
  crpsout_n[dd] <- mean(crps_norm(obs_today, loc, sc))
  crpsout_ens[dd] <- mean(crps_sample(obs_today, ensfc_today))
}

mean(crpsout_ens)
mean(crpsout_n)

# reasonably close to other model's results