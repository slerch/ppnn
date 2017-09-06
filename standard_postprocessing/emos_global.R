## benchmark EMOS model
rm(list=ls())

## load data
data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
library(ncdf4)
nc <- nc_open(paste0(data_dir, "data_interpolated.nc"))
fcdata <- ncvar_get(nc, "t2m_fc") 
## correct for forgetting to convert to degrees Celsius for years 2008-16
fcdata[fcdata > 50] <- fcdata[fcdata > 50] - 273.15
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
# version 0.9.3 or newer is required for exported parametric CRPS functions 
if(as.character(packageVersion("scoringRules")) < "0.9.3"){
  print("install scoringRules package version 0.9.3 or newer!")
}

# objective function for minimum CRPS estimation of parameters
# we fit an EMOS model N(mu, sigma), where
# the mean value mu = a + b * ensmean, and
# the variance sigma^2 = c + d * ensvar
# a,b,c,d are determined by minimizing the mean CRPS over a
# training set of past forecasts and observations
objective_fun <- function(par, ensmean, ensvar, obs){
  m <- cbind(1, ensmean) %*% par[1:2]
  ssq_tmp <- cbind(1, ensvar) %*% par[3:4]
  if(any(ssq_tmp < 0)){
    return(999999)
  } else{
    s <- sqrt(ssq_tmp)
    return(sum(crps_norm(y = obs, location = m, scale = s)))
  }
}

# gradient of objective function to use in optim()
#   including the gradient of the CRPS (as function of a,b,c,d)
#   allows faster numerical optmization as the gradient does
#   not need to be determined numerically
gradfun_wrapper <- function(par, obs, ensmean, ensvar){
  loc <- cbind(1, ensmean) %*% par[1:2]
  sc <- sqrt(cbind(1, ensvar) %*% par[3:4])
  dcrps_dtheta <- gradcrps_norm(y = obs, location = loc)
  out1 <- dcrps_dtheta[,1] %*% cbind(1, ensmean)
  out2 <- dcrps_dtheta[,2] %*% 
    cbind(1/(2*sqrt(par[3]+par[4]*ensvar)), 
          ensvar/(2*sqrt(par[3]+par[4]*ensvar)))
  return(as.numeric(cbind(out1,out2)))
}

# vector of validhours (00/12 UTC) to be able to select only forecast cases 
# with that specific valid hour for training
dates_validhour <- format(dates, '%H')

## main postprocessing function for a global EMOS 
# i.e., using all available stations for training. Argiments:
## input:
# vdate = validdate for which EMOS coefficients should be determined (should be from dates vector)
#         by the function
#         for 00 UTC forecasts, only past 00 UTC fc cases are used, similarly for 12 UTC
# train_length = length of rolling training period (in days)
## output: vector of coefficients a,b,c,d to be used at vdate
postproc_global <- function(vdate, train_length){
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
  
  # determine optimal EMOS coefficients a,b,c,d using minimum CRPS estimation
  optim_out <- optim(par = c(1,1,1,1), fn = objective_fun,
                     gr = gradfun_wrapper,
                     ensmean = ensfc_train_mean, ensvar = ensfc_train_var, 
                     obs = obs_train,
                     method = "BFGS")
  
  # check convergence of the numerical optimization
  if(optim_out$convergence != 0){
    message("numerical optimization did not converge")
  }

  # return optimal parameters
  return(optim_out$par)
}

## example for one date 
# tt <- as.POSIXct("2013-11-25 12:00", tz = "UTC", origin = "1970-01-01 00:00")
# m <- 25
# par_out <- postproc_global(vdate = tt, train_length = m)
# 
# ind_today <- which(dates == tt)
# ensfc_today <- fcdata[,,ind_today]
# ensfc_today_mean <- apply(ensfc_today, c(1), mean)
# ensfc_today_var <- apply(ensfc_today, c(1), var)
# 
# loc <- c(cbind(1, ensfc_today_mean) %*% par_out[1:2])
# sc <- sqrt(c(cbind(1, ensfc_today_var) %*% par_out[3:4]))
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

## example for 2008-2016, takes ~ 45 minutes on my laptop
# training length: 25 days
m <- 25
# dates for which EMOS models should be fit (here: only for 00 UTC valid times; all dates in 2008-2016)
dates_fit <- seq(as.POSIXct("2008-01-01 00:00", tz = "UTC", origin = "1970-01-01 00:00"), 
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
  par_out <- postproc_global(vdate = day, train_length = m)
  
  # extract ensemble forecasts and summary statistics for the day of interest ("day")
  ind_today <- which(dates == day)
  ensfc_today <- fcdata[,,ind_today]
  ensfc_today_mean <- apply(ensfc_today, c(1), mean)
  ensfc_today_var <- apply(ensfc_today, c(1), var)
  
  # use summary statistics and EMOS coefficients to determine mu and sigma
  # (= location and scale)
  # note that the link functions are the same as employed in the objective function
  loc <- c(cbind(1, ensfc_today_mean) %*% par_out[1:2])
  sc <- sqrt(c(cbind(1, ensfc_today_var) %*% par_out[3:4]))
  
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

summary(crpsout_n)
#   Min.  1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5874  0.8945  1.0198  1.0654  1.1827  4.4714 
summary(crpsout_ens)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5773  0.9726  1.1574  1.2175  1.3858  5.6086 

produce_pdfs <- FALSE
pdf_folder <- "/home/sebastian/Projects/PP_NN/code/standard_postprocessing/preliminary_results/"


## plot 1: time series of daily mean CRPS differences
if(produce_pdfs){
  pdf(paste0(pdf_folder, "crpsdiff_emos_global.pdf"), width = 6, height = 5, pointsize = 12)
}
plot(dates_fit, crpsout_ens - crpsout_n, type = "l")
abline(h = 0, lty = 2)
if(produce_pdfs){
  dev.off()
}

# months_fit <- format(dates_fit, "%Y%m")
# crpsout_ens_monthly <- NULL
# crpsout_n_monthly <- NULL
# for(i in 1:length(unique(months_fit))){
#   crpsout_ens_monthly[i] <- mean(crpsout_ens[which(months_fit == unique(months_fit)[i])])
#   crpsout_n_monthly[i] <- mean(crpsout_n[which(months_fit == unique(months_fit)[i])])
# }
# plot(crpsout_ens_monthly - crpsout_n_monthly, type = "l",
#      ylim = c(0,0.5))
# abline(v = 0:9*12, lty = 3)
# abline(h = 0, lty = 2)


## rolling mean
# based on http://www.cookbook-r.com/Manipulating_data/Calculating_a_moving_average/
movingAverage_ignoreNA <- function(x, n = 1) {
  before <- floor((n-1)/2)
  after <- ceiling((n-1)/2)
  
  # Track the sum and count of number of non-NA items
  s <- rep(0, length(x))
  count <- rep(0, length(x))
  
  # Add the centered data 
  new <- x
  # Add to count list wherever there isn't a 
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new
  
  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new <- c(rep(NA, i), x[1:(length(x)-i)])
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    i <- i+1
  }
  
  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new   <- c(x[(i+1):length(x)], rep(NA, i))
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    i <- i+1
  }
  
  # return sum divided by count
  s/count
}


# plot 2: rolling 50-day mean (mean) CRPS values and differences (smoothed version of plot 1)
k <- 50
if(produce_pdfs){
  pdf(paste0(pdf_folder, "crps_rollingmean_emos_global.pdf"), width = 12, height = 5, pointsize = 12)
}
par(mfrow = c(1,2))
plot(dates_fit, movingAverage_ignoreNA(crpsout_ens, k), type = "l", ylim = c(0.75,1.75),
     main = "50 day rolling mean", ylab = "mean CRPS")
lines(dates_fit, movingAverage_ignoreNA(crpsout_n, k), col = "blue")
legend("topright", legend = c("Ensemble", "EMOS global"), lty = c(1,1), col = c("black", "blue"), bty = "n")

plot(dates_fit, movingAverage_ignoreNA(crpsout_ens, k)-movingAverage_ignoreNA(crpsout_n, k), type = "l", 
     main = "50 day rolling mean", ylab = "mean CRPS difference: Ens - EMOS")
abline(h = 0, lty = 2)
if(produce_pdfs){
  dev.off()
}

# plot 3: separate forecast cases by month and average over all years
months_fit <- format(dates_fit, "%m")
crpsout_ens_monthly <- NULL
crpsout_n_monthly <- NULL
for(i in 1:length(unique(months_fit))){
  crpsout_ens_monthly[i] <- mean(crpsout_ens[which(months_fit == unique(months_fit)[i])])
  crpsout_n_monthly[i] <- mean(crpsout_n[which(months_fit == unique(months_fit)[i])])
}
# plot(unique(format(dates_fit, "%m")), crpsout_ens_monthly - crpsout_n_monthly, type = "l",
#      ylim = c(0,0.3))
# abline(h = 0, lty = 2)

if(produce_pdfs){
  pdf(paste0(pdf_folder, "crps_monthlymean_emos_global.pdf"), width = 6, height = 5, pointsize = 12)
}
plot(unique(format(dates_fit, "%m")), crpsout_ens_monthly, type = "o", ylim = range(c(crpsout_ens_monthly, crpsout_n_monthly)),
     main = "monthly mean CRPS", ylab = "mean CRPS", xlab = "month")
lines(unique(format(dates_fit, "%m")), crpsout_n_monthly, col = "blue", type = "o")
legend("top", legend = c("Ensemble", "EMOS global"), lty = c(1,1), col = c("black", "blue"), bty = "n")
if(produce_pdfs){
  dev.off()
}