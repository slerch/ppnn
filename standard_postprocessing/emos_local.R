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

dates_validhour <- format(dates, '%H')

# local EMOS 
# date = validdate for which EMOS coefficients should be determined (should be from dates vector)
#         for 00 UTC forecasts, only past 00 UTC fc cases are used, similarly for 12 UTC
postproc_local <- function(vdate, train_length){
  
  # matrix to save output to
  par_return <- matrix(NA, nrow = length(stations), ncol = 4)
  
  # determine training set
  vhour <- format(vdate, "%H")
  
  # extract start and end date of rolling training period
  enddate <- vdate - 2*24*3600
  startdate <- enddate - (train_length-1)*24*3600
  if(startdate < dates[1]){
    stop("start date of training period before first available date, choose shorter train_length")
  }
  
  ind_training <- which(dates_validhour == vhour & dates >= startdate & dates <= enddate)
  
  # loop over stations
  for(this_station in stations){
    
    ind_this_station <- which(stations == this_station)
    
    # extract relevant subset of fcdata
    ensfc_train <- fcdata[ind_this_station,,ind_training]  
    ensfc_train_mean <- apply(ensfc_train, c(2), mean)
    ensfc_train_var <- apply(ensfc_train, c(2), var)
    
    # extract relevant subset of obsdata
    obs_train <- obsdata[ind_this_station,ind_training]
    
    # find and omit NA values
    if(anyNA(obs_train)){
      ind_notNA <- which(!is.na(obs_train))
      obs_train <- obs_train[ind_notNA]
      ensfc_train_mean <- ensfc_train_mean[ind_notNA]
      ensfc_train_var <- ensfc_train_var[ind_notNA]
    }
    
    # account for cases with only missing observations
    if(length(obs_train) == 0){
      par_return[ind_this_station,] <- rep(NA, 4)
    } else{
      optim_out <- optim(par = c(1,1,1,1), fn = objective_fun,
                         gr = gradfun_wrapper,
                         ensmean = ensfc_train_mean, ensvar = ensfc_train_var, 
                         obs = obs_train,
                         method = "BFGS")
      
      par_return[ind_this_station,] <- optim_out$par
    }
    
  }
  
  # return optimal parameters
  return(par_return)
}

## example for one date 
# tt <- as.POSIXct("2013-11-25 12:00", tz = "UTC", origin = "1970-01-01 00:00")
# m <- 125
# par_out <- postproc_local(vdate = tt, train_length = m)
# 
# ind_today <- which(dates == tt)
# ensfc_today <- fcdata[,,ind_today]
# ensfc_today_mean <- apply(ensfc_today, c(1), mean)
# ensfc_today_var <- apply(ensfc_today, c(1), var)
# 
# loc <- NULL
# sc <- NULL
# for(i in 1:nrow(par_out)){
#   loc[i] <- c(c(1, ensfc_today_mean[i]) %*% par_out[i,1:2])
#   # correct for possible negative scale values (only a preliminary fix here; should be done better (e.g. using L-BFGS-B))
#   tmp <- sqrt(c(c(1, ensfc_today_var[i]) %*% par_out[i,3:4]))
#   if(!is.finite(tmp)){
#     tmp <- sqrt(abs(c(c(1, ensfc_today_var[i]) %*% par_out[i,3:4])))
#   }
#   sc[i] <- tmp
# }
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

## --- apply local EMOS model to 2016 for different train_length values
##      takes ~ 12 hours on my laptop

# m_values <- c(50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 700)
# dates_fit <- seq(as.POSIXct("2016-01-01 00:00", tz = "UTC", origin = "1970-01-01 00:00"), 
#                  as.POSIXct("2016-12-31 00:00", tz = "UTC", origin = "1970-01-01 00:00"), 
#                  by="+1 day")
# crpsout_n <- matrix(nrow = length(dates_fit), ncol = length(m_values))
# 
# for(dd in 1:length(dates_fit)){
#   day <- dates_fit[dd]
#   if(format(dates_fit[dd], "%d") == "01"){
#     format(dates_fit[dd], "%B %Y")
#     cat("Starting at", paste(Sys.time()), ":", format(dates_fit[dd], "%B %Y"), "\n"); flush(stdout())
#   }
#   
#   ind_today <- which(dates == day)
#   ensfc_today <- fcdata[,,ind_today]
#   ensfc_today_mean <- apply(ensfc_today, c(1), mean)
#   ensfc_today_var <- apply(ensfc_today, c(1), var)
#   
#   for(mm in m_values){
#     par_out <- postproc_local(vdate = day, train_length = mm)
#     
#     loc <- NULL
#     sc <- NULL
#     for(i in 1:nrow(par_out)){
#       loc[i] <- c(c(1, ensfc_today_mean[i]) %*% par_out[i,1:2])
#       # correct for possible negative scale values (only a preliminary fix here; should be done better (e.g. using L-BFGS-B))
#       tmp <- sqrt(c(c(1, ensfc_today_var[i]) %*% par_out[i,3:4]))
#       if(!is.finite(tmp)){
#         tmp <- sqrt(abs(c(c(1, ensfc_today_var[i]) %*% par_out[i,3:4])))
#       }
#       sc[i] <- tmp
#     }
#     
#     obs_today <- obsdata[,ind_today]
#     
#     if(anyNA(c(obs_today, loc, sc))){
#       ind_notNA_today_obs <- which(!is.na(obs_today))
#       ind_notNA_today_loc <- which(!is.na(loc))
#       ind_notNA_today_sc <- which(!is.na(sc))
#       ind_notNA_today <- intersect(intersect(ind_notNA_today_obs, ind_notNA_today_loc), ind_notNA_today_sc)
#       obs_today <- obs_today[ind_notNA_today]
#       loc <- loc[ind_notNA_today]
#       sc <- sc[ind_notNA_today]
#     }
#     
#     crpsout_n[dd, which(m_values == mm)] <- mean(crps_norm(obs_today, loc, sc))
#   }
# }
# 
# summary(crpsout_n)
# 
# plot(m_values, apply(crpsout_n, 2, mean))


## --- entire period

## example for 2008-2016, takes ~ 10 hours on my laptop
m <- 360
dates_fit <- seq(as.POSIXct("2008-01-01 00:00", tz = "UTC", origin = "1970-01-01 00:00"), 
                 as.POSIXct("2016-12-31 00:00", tz = "UTC", origin = "1970-01-01 00:00"), 
                 by="+1 day")
crpsout_n <- numeric(length(dates_fit))
crpsout_ens <- numeric(length(dates_fit))

for(dd in 1:length(dates_fit)){
  day <- dates_fit[dd]
  if(format(dates_fit[dd], "%d") == "01"){
    format(dates_fit[dd], "%B %Y")
    cat("Starting at", paste(Sys.time()), ":", format(dates_fit[dd], "%B %Y"), "\n"); flush(stdout())
  }
  par_out <- postproc_local(vdate = day, train_length = m)
  
  ind_today <- which(dates == day)
  ensfc_today <- fcdata[,,ind_today]
  ensfc_today_mean <- apply(ensfc_today, c(1), mean)
  ensfc_today_var <- apply(ensfc_today, c(1), var)
  
  loc <- NULL
  sc <- NULL
  for(i in 1:nrow(ensfc_today)){
    loc[i] <- c(c(1, ensfc_today_mean[i]) %*% par_out[i,1:2])
    # correct for possible negative scale values (only a preliminary fix here; should be done better (e.g. using L-BFGS-B))
    tmp <- sqrt(c(c(1, ensfc_today_var[i]) %*% par_out[i,3:4]))
    if(!is.finite(tmp)){
      tmp <- sqrt(abs(c(c(1, ensfc_today_var[i]) %*% par_out[i,3:4])))
    }
    sc[i] <- tmp
  }
  
  obs_today <- obsdata[,ind_today]
  
  if(anyNA(obs_today)){
    ind_notNA_today <- which(!is.na(obs_today))
    obs_today <- obs_today[ind_notNA_today]
    loc <- loc[ind_notNA_today]
    sc <- sc[ind_notNA_today]
    ensfc_today <- ensfc_today[ind_notNA_today,]
  }
  
  crpsout_n[dd] <- mean(crps_norm(obs_today, loc, sc))
  crpsout_ens[dd] <- mean(crps_sample(obs_today, ensfc_today))
}

### results local EMOS model
summary(crpsout_n)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.4960  0.7673  0.9008  0.9578  1.0756  5.2003      73 
### compare to globale results:
#   Min.  1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5874  0.8945  1.0198  1.0654  1.1827  4.4714 
### results ensemble
summary(crpsout_ens)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.5773  0.9726  1.1574  1.2175  1.3858  5.6086 

## save results
# save(crpsout_n, crpsout_ens,
#      file = "/home/sebastian/Projects/PP_NN/code/standard_postprocessing/preliminary_results/scores_local.Rdata")


## plots

produce_pdfs <- TRUE
pdf_folder <- "/home/sebastian/Projects/PP_NN/code/standard_postprocessing/preliminary_results/"


## plot 1: time series of daily mean CRPS differences
if(produce_pdfs){
  pdf(paste0(pdf_folder, "crpsdiff_emos_local.pdf"), width = 6, height = 5, pointsize = 12)
}
plot(dates_fit, crpsout_ens - crpsout_n, type = "l")
abline(h = 0, lty = 2)
if(produce_pdfs){
  dev.off()
}



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
  pdf(paste0(pdf_folder, "crps_rollingmean_emos_local.pdf"), width = 12, height = 5, pointsize = 12)
}
par(mfrow = c(1,2))
plot(dates_fit, movingAverage_ignoreNA(crpsout_ens, k), type = "l", ylim = c(0.75,1.75),
     main = "50 day rolling mean", ylab = "mean CRPS")
lines(dates_fit, movingAverage_ignoreNA(crpsout_n, k), col = "blue")
legend("topright", legend = c("Ensemble", "EMOS local"), lty = c(1,1), col = c("black", "blue"), bty = "n")

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
  crpsout_n_monthly[i] <- mean(crpsout_n[which(months_fit == unique(months_fit)[i])], na.rm = TRUE)
}
# plot(unique(format(dates_fit, "%m")), crpsout_ens_monthly - crpsout_n_monthly, type = "l",
#      ylim = c(0,0.3))
# abline(h = 0, lty = 2)

if(produce_pdfs){
  pdf(paste0(pdf_folder, "crps_monthlymean_emos_local.pdf"), width = 6, height = 5, pointsize = 12)
}
plot(unique(format(dates_fit, "%m")), crpsout_ens_monthly, type = "o", ylim = range(c(crpsout_ens_monthly, crpsout_n_monthly)),
     main = "monthly mean CRPS", ylab = "mean CRPS", xlab = "month")
lines(unique(format(dates_fit, "%m")), crpsout_n_monthly, col = "blue", type = "o")
legend("top", legend = c("Ensemble", "EMOS local"), lty = c(1,1), col = c("black", "blue"), bty = "n")
if(produce_pdfs){
  dev.off()
}
