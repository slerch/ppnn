## run selected benchmark models
## for documented code see model run code in ../standard_postprocessing/ and ../qrf_postprocessing/
## see note.md

## training 2007-2015
## EMOS local, fixed window

## chosen tuning parameters:
##    none


rm(list=ls())

data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
load(paste0(data_dir, "data_all.Rdata"))
data <- data[, -which(!(names(data) %in% c("obs", "date", "station", "t2m_mean", "t2m_var")))]

library(scoringRules)
library(lubridate)

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

stations_list <- unique(data$station)

train_end <- as.Date("2016-01-01 00:00 UTC") - days(2)
train_start <- data$date[1]

data_train_all <- subset(data, date >= train_start & date <= train_end)

eval_start <- as.Date("2016-01-01 00:00 UTC")
eval_end <- as.Date("2016-12-31 00:00 UTC")
eval_dates <- seq(eval_start, eval_end, by = "1 day")

data_eval_all <- subset(data, date >= eval_start & date <= eval_end)

out_loc <- NULL
out_sc <- NULL


t1 <- Sys.time()

out_loc <- rep(NA, nrow(data_eval_all))
out_sc <- rep(NA, nrow(data_eval_all))

par_out <- matrix(NA, ncol = 4, nrow = length(stations_list))

for(this_station in stations_list){
  if(which(stations_list == this_station) %% 10 == 0){print(which(stations_list == this_station))}
  
  data_train <- subset(data_train_all, station == this_station)
  
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
  
  data_eval <- subset(data_eval_all, station == this_station)
  
  loc_st <- c(cbind(1, data_eval$t2m_mean) %*% optim_out$par[1:2])
  scsquared_tmp <- c(cbind(1, data_eval$t2m_var) %*% optim_out$par[3:4])
  sc_st <- sqrt(abs(scsquared_tmp))
  
  # save parameters
  ind_this_station <- which(data_eval_all$station == this_station)
  
  out_loc[ind_this_station] <- loc_st 
  out_sc[ind_this_station] <- sc_st
}

t2 <- Sys.time()
t2-t1
# run time ~ 43 seconds

df_out <- data.frame(date = data_eval_all$date,
                     station_id = data_eval_all$station,
                     mean = out_loc,
                     std = out_sc)

summary(crps_norm(data_eval_all$obs, mean = out_loc, sd = out_sc))

emos_local_fixed_alldata <- df_out
save(emos_local_fixed_alldata,
     file = "/home/sebastian/Projects/PP_NN/code/results/Rdata_files/emos_local_fixed_2007-15.Rdata")
