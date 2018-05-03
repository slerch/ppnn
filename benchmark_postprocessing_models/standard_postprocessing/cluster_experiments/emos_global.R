## implementation of standard EMOS model with rolling training window
## here: based on prepared Rdata file

## --------------- data preparation --------------- ##

rm(list=ls())

data_dir <- "..."
load(paste0(data_dir, "data_all.Rdata"))

# renove all covariates (not used in standard EMOS)
data <- data[, -which(!(names(data) %in% c("obs", "date", "station", "t2m_mean", "t2m_var")))]
head(data)

library(scoringRules, lib = "...")
library(lubridate, lib = "...")

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


## --------------- post-processing function --------------- ##

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
  train_end <- vdate - days(2)
  train_start <- train_end - days(train_length - 1)
  
  data_train <- subset(data, date >= train_start & date <= train_end)
  
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


## --------------- wrapper function for execution on cluster --------------- ##

eval_start <- as.Date("2016-01-01 00:00 UTC")
eval_end <- as.Date("2016-12-31 00:00 UTC")
eval_dates <- seq(eval_start, eval_end, by = "1 day")
data_eval_all <- subset(data, date >= eval_start & date <= eval_end)

run_pp <- function(m){
  
  Routname1 <- paste0("emos_global_m", m, ".Rout")
  Routname <- paste0("...", Routname1)
  sink(Routname)
  
  crps_pp <- NULL
  
  for(day_id in 1:length(eval_dates)){
    
    today <- eval_dates[day_id]
    
    # progress indicator
    if(day(as.Date(today)) == 1){
      cat("Starting at", paste(Sys.time()), ":", as.character(today), "\n"); flush(stdout())
    }
    
    # post-processing
    par_out <- postproc_global(vdate = today, train_length = m)
    
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
    crps_today <- crps_norm(y = data_eval$obs, mean = loc, sd = sc)
    crps_pp[which(data_eval_all$date == today)] <- crps_today
  }
  
  savename <- paste0("emos_global_m", m, ".Rdata")
  fname <- paste0("...", savename)
  save(crps_pp, file = fname)
  
  sink()
}

m_val <- seq(5, 365, by = 5) # 73 tasks

ID <- as.numeric(as.character(Sys.getenv(c("SGE_TASK_ID"))))

run_pp(m = m_val[ID])