## implementation of local EMOS model with boosting 
## here: code to run experiments on cluster to test tuning parameters
## here: based on prepared Rdata file

## --------------- data preparation --------------- ##

rm(list=ls())

data_dir <- "..."
load(paste0(data_dir, "data_all.Rdata"))

data$sm_mean <- NULL
data$sm_var <- NULL

library(scoringRules, lib = "...")
library(lubridate, lib = "...")
library(crch, lib = "...")

train_end <- as.Date("2016-01-01 00:00 UTC") - days(2)
train_start <- data$date[1]

data_train_all <- subset(data, date >= train_start & date <= train_end)

eval_start <- as.Date("2016-01-01 00:00 UTC")
eval_end <- as.Date("2016-12-31 00:00 UTC")
eval_dates <- seq(eval_start, eval_end, by = "1 day")

data_eval_all <- subset(data, date >= eval_start & date <= eval_end)

run_boosting <- function(maxit, nu, mstop, linkfc){
  
  Routname <- paste0("...",
                     "alldata_boosting_local_maxit", maxit,
                     "_nu", nu,
                     "_mstop-", mstop,
                     "_link-", linkfc,
                     ".Rout")
  sink(Routname)
  
  loc_save <- numeric(nrow(data_eval_all))
  sc_save <- numeric(nrow(data_eval_all))
  
  stations_list <- unique(data$station)
  
  for(this_station in stations_list){
    # progress indicator
    progind <- which(stations_list == this_station)
    if(progind %% 10 == 0){
      cat(progind, "of", length(stations_list), "started at", paste(Sys.time()), "\n")
    }
    
    # data_train <- subset(data, date >= train_start & date <= train_end & station == this_station)
    data_train <- subset(data_train_all, station == this_station)
    
    # remove incomplete cases (= NA obs or fc)
    data_train <- data_train[complete.cases(data_train), ]
    if(nrow(data_train) < 10){
      next
    }
    
    crch_model <- crch(obs ~ .|.,
                       data = data_train[,-which(names(data) %in% c("date", "station"))],
                       dist = "gaussian",
                       link.scale = linkfc,
                       nu = nu,
                       method = "boosting",
                       maxit = maxit,
                       mstop = mstop)
    
    data_eval <- subset(data_eval_all, station == this_station)
    
    loc <- as.numeric(predict(crch_model, data_eval, type = "location"))
    sc <- as.numeric(predict(crch_model, data_eval, type = "scale"))
    
    # save parameters
    ind_this_station <- which(data_eval_all$station == this_station)
    
    loc_save[ind_this_station] <- loc 
    sc_save[ind_this_station] <- sc
  }
  
  crps_boosting <- crps_norm(y = data_eval_all$obs, mean = loc_save, sd = sc_save)
  
  savename <- paste0("alldata_boosting_local_maxit", maxit,
                     "_nu", nu,
                     "_mstop-", mstop,
                     "_link-", linkfc,
                     ".Rdata")
  fname <- paste0("...", savename)
  save(crps_boosting, file = fname)
  
  print(mean(crps_boosting, na.rm = TRUE))
  
  sink()
}

## tuning parameters for boosting:
# maxit = number of boosting iterations, default = 100
# nu = boosting step size
# mstop = c("max", "aic", "bic", "cv")
#     nfolds = number of folds if mstop = cv
# link.scale = c("log", "identity", "quadratic")

maxit_try <- c(250, 500, 1000, 2000, 5000)
nu_try <- c(0.1, 0.025, 0.05, 0.2, 0.4)
mstop_try <- c("max", "aic", "bic", "cv")
link_try <-  c("log")

pars <- expand.grid(maxit_try, nu_try, mstop_try, link_try)
nrow(pars) # 100 tasks

ID <- as.numeric(as.character(Sys.getenv(c("SGE_TASK_ID"))))

run_boosting(maxit = pars[ID,1], 
             nu = pars[ID,2], 
             mstop = as.character(pars[ID,3]), 
             linkfc = as.character(pars[ID,4]))