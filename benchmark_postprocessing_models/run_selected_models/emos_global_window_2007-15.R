## run selected benchmark models
## for documented code see model run code in ../standard_postprocessing/ and ../qrf_postprocessing/
## see note.md

## training 2007-2015
## EMOS global, "climatological" window: for day d0, use [d0 - m days, d0 + m days] in all previous years

## chosen tuning parameters:
##    m = 20 days


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


postproc_global <- function(vdate, window_width_half){
  
  # determine training set
  # deal with February 29 existing only all 4 years (set to Feb 28)
  if(month(vdate) == 2){
    if(day(vdate) == 29){
      vdate <- vdate - days(1)
    }
  }
  
  vdate_prevyears <- vdate - years(1:10)
  if(any(vdate_prevyears < data$date[1])){
    vdate_prevyears[which(vdate_prevyears < data$date[1])] <- NULL
  }
  
  train_dates <- vdate_prevyears
  for(yy in 1:length(vdate_prevyears)){
    today <- vdate_prevyears[yy]
    train_dates <- c(train_dates, today - days(1:window_width_half))
    train_dates <- c(train_dates, today + days(1:window_width_half))
  }
  
  data_train <- data[which(is.element(as.Date(data$date), as.Date(train_dates))), ]
  
  # remove incomplete cases (= NA obs or fc)
  data_train <- data_train[complete.cases(data_train), ]
  
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
  
  # return optimal parameters
  return(optim_out$par)
}

eval_start <- as.Date("2016-01-01 00:00 UTC")
eval_end <- as.Date("2016-12-31 00:00 UTC")
eval_dates <- seq(eval_start, eval_end, by = "1 day")

m <- 20

out_loc <- NULL
out_sc <- NULL

data_eval_all <- subset(data, date >= eval_start & date <= eval_end)

t1 <- Sys.time()

for(day_id in 1:length(eval_dates)){
  
  today <- eval_dates[day_id]
  
  # progress indicator
  if(day(as.Date(today)) == 1){
    cat("Starting at", paste(Sys.time()), ":", as.character(today), "\n"); flush(stdout())
  }
  
  # post-processing
  par_out <- postproc_global(vdate = today, window_width_half = m)
  
  # out of sample distribution parameters for today
  data_eval <- subset(data, date == today)
  loc <- c(cbind(1, data_eval$t2m_mean) %*% par_out[1:2])
  scsquared_tmp <- c(cbind(1, data_eval$t2m_var) %*% par_out[3:4])
  if(any(scsquared_tmp <= 0)){
    print("negative scale, taking absolute value")
    sc <- sqrt(abs(scsquared_tmp))
  } else{
    sc <- sqrt(scsquared_tmp)
  }
  
  # CRPS computation
  out_loc <- c(out_loc, loc)
  out_sc <- c(out_sc, sc)
}

t2 <- Sys.time()
t2-t1
# run time ~ 12 minutes

df_out <- data.frame(date = data_eval_all$date,
                     station_id = data_eval_all$station,
                     mean = out_loc,
                     std = out_sc)

summary(crps_norm(data_eval_all$obs, out_loc, out_sc))

emos_global_window_alldata <- df_out
save(emos_global_window_alldata,
     file = "/home/sebastian/Projects/PP_NN/code/results/Rdata_files/emos_global_window_2007-15.Rdata")
