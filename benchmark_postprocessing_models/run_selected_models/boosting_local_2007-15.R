## run selected benchmark models
## for documented code see model run code in ../standard_postprocessing/ and ../qrf_postprocessing/
## see note.md

## training 2007-15
## EMOS+Boosting local

## chosen tuning parameters (same as 2015 training):
##    maxit = 1000, nu = 0.05, mstop = "aic"
##    (results are robust to wide variety of tuning parameter choices)

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
train_start <- data$date[1]

data_train_all <- subset(data, date >= train_start & date <= train_end)

eval_start <- as.Date("2016-01-01 00:00 UTC")
eval_end <- as.Date("2016-12-31 00:00 UTC")
eval_dates <- seq(eval_start, eval_end, by = "1 day")

data_eval_all <- subset(data, date >= eval_start & date <= eval_end)

out_loc <- rep(NA, nrow(data_eval_all))
out_sc <- rep(NA, nrow(data_eval_all))

stations_list <- unique(data$station)

t1 <- Sys.time()

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
  
  ## NOTE: boosting is only implemented for link.scale = "log", otherwise cryptic error message 
  crch_model <- crch(obs ~ .|.,
                     data = data_train[,-which(names(data) %in% c("date", "station"))],
                     dist = "gaussian",
                     link.scale = "log",
                     method = "boosting",
                     maxit = 1000,
                     nu = 0.05,
                     mstop = "aic")
  
  data_eval <- subset(data_eval_all, station == this_station)
  
  loc <- as.numeric(predict(crch_model, data_eval, type = "location"))
  sc <- as.numeric(predict(crch_model, data_eval, type = "scale"))
  
  # save parameters
  ind_this_station <- which(data_eval_all$station == this_station)
  
  out_loc[ind_this_station] <- loc 
  out_sc[ind_this_station] <- sc
}

t2 <- Sys.time()
t2-t1
# run time ~ 48 minutes
# not directly comparable to std EMOS are due to implementation in crch vs direct optimization implementation and ML vs minCRPS
# fair comparions would be to use crch implementation for std EMOS model as well

df_out <- data.frame(date = data_eval_all$date,
                     station_id = data_eval_all$station,
                     mean = out_loc,
                     std = out_sc)

summary(crps_norm(data_eval_all$obs, mean = out_loc, sd = out_sc))

boosting_local_alldata <- df_out
save(boosting_local_alldata,
     file = "/home/sebastian/Projects/PP_NN/code/results/Rdata_files/boosting_local_2007-15.Rdata")
