## implementation of local EMOS model with boosting 
## training on data from 2015 only
## here: results for some selected tuning parameters
## here: based on prepared Rdata file

## --------------- data preparation --------------- ##

rm(list=ls())

data_dir <- "/media/sebastian/Elements/Postproc_NN/data/"
load(paste0(data_dir, "data_all.Rdata"))

# remove sm data (missing values)
data$sm_mean <- NULL
data$sm_var <- NULL
head(data)

library(scoringRules)
library(lubridate)
library(crch)

train_end <- as.Date("2016-01-01 00:00 UTC") - days(2)
train_start <- as.Date("2015-01-01 00:00 UTC") 

data_train_all <- subset(data, date >= train_start & date <= train_end)

eval_start <- as.Date("2016-01-01 00:00 UTC")
eval_end <- as.Date("2016-12-31 00:00 UTC")
eval_dates <- seq(eval_start, eval_end, by = "1 day")

data_eval_all <- subset(data, date >= eval_start & date <= eval_end)

## -------------- model estimation --------------  ##

## standard model without any boosting

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
                     link.scale = "log",
                     method = "boosting",
                     maxit = 100,
                     mstop = "aic")
  
  data_eval <- subset(data_eval_all, station == this_station)
  
  loc <- as.numeric(predict(crch_model, data_eval, type = "location"))
  sc <- as.numeric(predict(crch_model, data_eval, type = "scale"))
  
  # save parameters
  ind_this_station <- which(data_eval_all$station == this_station)
  
  loc_save[ind_this_station] <- loc 
  sc_save[ind_this_station] <- sc
}

crps_boost <- crps_norm(y = data_eval_all$obs, mean = loc_save, sd = sc_save)
summary(crps_boost)

