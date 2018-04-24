## implementation of standard EMOS model with rolling training window
## here: based on prepared Rdata file
## test version with non-rolling training window (all of year 2015)

## --------------- data preparation --------------- ##

rm(list=ls())

data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
load(paste0(data_dir, "data_all.Rdata"))

# renove all covariates (not used in standard EMOS)
# data[, !(names(data) %in% c("obs", "date", "station", "t2m_mean", "t2m_var"))] <- NULL # does not work on cluster?
data <- data[, -which(!(names(data) %in% c("obs", "date", "station", "t2m_mean", "t2m_var")))]
head(data)

library(scoringRules)
library(lubridate)

## --------------- Helper functions for optimization --------------- ##

# objective function for minimum CRPS estimation of parameters
#   we fit an EMOS model N(mu, sigma), where
#   the mean value mu = a + b * ensmean, and
#   the variance sigma^2 = c + d * ensvar
#   a,b,c,d are determined by minimizing the mean CRPS over a
#   training set of past forecasts and observations
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
  dcrps_dtheta <- gradcrps_norm(y = obs, location = loc, scale = sc) 
  out1 <- dcrps_dtheta[,1] %*% cbind(1, ensmean)
  out2 <- dcrps_dtheta[,2] %*% 
    cbind(1/(2*sqrt(par[3]+par[4]*ensvar)), 
          ensvar/(2*sqrt(par[3]+par[4]*ensvar)))
  return(as.numeric(cbind(out1,out2)))
}

## --------------- post-processing  --------------- ##
## post-processing

stations_list <- unique(data$station)

# determine training set
train_end <- as.Date("2016-01-01 00:00 UTC") - days(2)
train_start <- data$date[1]

data_train_dates <- subset(data, date >= train_start & date <= train_end)

par_out <- matrix(NA, ncol = 4, nrow = length(stations_list))

# loop over stations
for(this_station in stations_list){
  if(which(stations_list == this_station) %% 10 == 0){print(which(stations_list == this_station))}
  
  # data_train <- subset(data, date >= train_start & date <= train_end & station == this_station)
  data_train <- subset(data_train_dates, station == this_station)
  
  # remove incomplete cases (= NA obs or fc)
  data_train <- data_train[complete.cases(data_train), ]
  
  # skip station if there are too few forecast cases
  if(nrow(data_train) < 10){
    next
  }
  
  # determine optimal EMOS coefficients a,b,c,d using minimum CRPS estimation
  optim_out <- optim(par = c(1,1,1,1), 
                     fn = objective_fun,
                     gr = gradfun_wrapper,
                     ensmean = data_train$t2m_mean, 
                     ensvar = data_train$t2m_var, 
                     obs = data_train$obs,
                     method = "BFGS")
  
  # check convergence of the numerical optimization
  if(optim_out$convergence != 0){
    message("numerical optimization did not converge")
  }
  
  par_out[which(stations_list == this_station), ] <- optim_out$par
}


# 2016
eval_start <- as.Date("2016-01-01 00:00 UTC")
eval_end <- as.Date("2016-12-31 00:00 UTC")
eval_dates <- seq(eval_start, eval_end, by = "1 day")

crps_pp <- NULL

data_eval_all <- subset(data, date >= eval_start & date <= eval_end)

for(day_id in 1:length(eval_dates)){
  
  today <- eval_dates[day_id]
  
  # progress indicator
  if(day(as.Date(today)) == 1){
    cat("Starting at", paste(Sys.time()), ":", as.character(today), "\n"); flush(stdout())
  }
  
  # out of sample distribution parameters for today
  loc <- rep(NA, length(stations_list))
  sc <- rep(NA, length(stations_list))
  for(this_station in stations_list){
    ind_st <- which(stations_list == this_station)
    # print(ind_st)
    data_eval <- subset(data_eval_all, date == today & station == this_station)
    if(!is.finite(data_eval$obs)){next}
    loc_st <- c(cbind(1, data_eval$t2m_mean) %*% par_out[ind_st,1:2])
    scsquared_tmp <- c(cbind(1, data_eval$t2m_var) %*% par_out[ind_st,3:4])
    if(is.na(scsquared_tmp)){
      next
    }
    if(scsquared_tmp <= 0){
      print("negative scale, taking absolute value")
      sc_st <- sqrt(abs(scsquared_tmp))
    } else{
      sc_st <- sqrt(scsquared_tmp)
    }
    loc[ind_st] <- loc_st
    sc[ind_st] <- sc_st
  }
  
  # CRPS
  data_eval_today <- subset(data_eval_all, date == today)
  crps_today <- crps_norm(y = data_eval_today$obs, mean = loc, sd = sc)
  crps_pp[which(data_eval_all$date == today)] <- crps_today
}

summary(crps_pp)

